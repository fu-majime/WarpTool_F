//! AviUtl2 script module for contour-based puppet meshes.
//!
//! The module accepts the data pointer returned by `obj.getpixeldata()` and
//! returns ordinary Lua arrays. No allocator-owned Rust memory crosses the
//! plugin boundary.

use std::{
    ffi::CString,
    panic::{AssertUnwindSafe, catch_unwind},
    ptr::NonNull,
};

use aviutl2::{
    AnyResult,
    anyhow::{self, Context as _},
    module::{ModuleFunction, ParamType, ScriptModuleCallHandle, ScriptModuleFunctions},
    sys::module2::SCRIPT_MODULE_PARAM,
};

mod contour;
mod delaunay;
mod dilate;
mod poisson;
mod simplify;

const RGBA_CHANNELS: usize = 4;
const MIN_DENSITY: i32 = 5;
const MAX_DENSITY: i32 = 200;

#[aviutl2::plugin(ScriptModule)]
struct PuppetMeshModule;

impl aviutl2::module::ScriptModule for PuppetMeshModule {
    fn new(_info: aviutl2::AviUtl2Info) -> AnyResult<Self> {
        Ok(Self)
    }

    fn plugin_info(&self) -> aviutl2::module::ScriptModuleTable {
        aviutl2::module::ScriptModuleTable {
            information: format!(
                "Puppet mesh generator for AviUtl2 / v{} / Fu-Majime",
                env!("CARGO_PKG_VERSION")
            ),
            functions: Self::functions(),
        }
    }
}

impl ScriptModuleFunctions for PuppetMeshModule {
    fn functions() -> Vec<ModuleFunction> {
        vec![ModuleFunction {
            name: "generate".to_owned(),
            func: generate_callback,
        }]
    }
}

aviutl2::register_script_module!(PuppetMeshModule);

/// AviUtl2 exposes the value returned by `obj.getpixeldata()` as Lua
/// `Userdata`. aviutl2-rs 0.42's typed pointer conversion only accepts
/// `LightUserdata`, so this one callback deliberately reads the pointer through
/// the SDK's raw `get_param_data` entry point. All other arguments and results
/// still use the checked high-level wrapper.
extern "C" fn generate_callback(raw: *mut SCRIPT_MODULE_PARAM) {
    let result = catch_unwind(AssertUnwindSafe(|| generate_callback_inner(raw)));

    match result {
        Ok(Ok(())) => {}
        Ok(Err(error)) => set_callback_error(raw, &error.to_string()),
        Err(payload) => {
            let message = payload
                .downcast_ref::<&str>()
                .copied()
                .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
                .unwrap_or("unknown panic");
            set_callback_error(raw, &format!("internal panic: {message}"));
        }
    }
}

fn generate_callback_inner(raw: *mut SCRIPT_MODULE_PARAM) -> AnyResult<()> {
    anyhow::ensure!(!raw.is_null(), "AviUtl2 passed a null module-call handle");

    // SAFETY: AviUtl2 owns `raw` and guarantees that it remains valid for the
    // duration of this synchronous module callback.
    let mut params = unsafe { ScriptModuleCallHandle::from_raw(raw) };
    let data_type = params.get_param_type(0);
    anyhow::ensure!(
        matches!(
            data_type,
            Some(ParamType::Userdata | ParamType::LightUserdata)
        ),
        "parameter #0 must be pixel data (Userdata), got {data_type:?}"
    );

    // SAFETY: `raw` was checked above and points to a valid SDK parameter
    // table. Unlike aviutl2-rs 0.42's wrapper, the SDK function supports the
    // full Userdata value returned by `obj.getpixeldata()`.
    let data = unsafe { ((*raw).get_param_data)(0) };
    let data = NonNull::new(data.cast::<u8>()).context("pixel data pointer is null")?;

    let width = params.get_param_int(1).context("invalid width parameter")?;
    let height = params
        .get_param_int(2)
        .context("invalid height parameter")?;
    let threshold = params
        .get_param_int(3)
        .context("invalid threshold parameter")?;
    let density = params
        .get_param_int(4)
        .context("invalid density parameter")?;
    let border_px = params
        .get_param_int(5)
        .context("invalid border parameter")?;

    let (vertices, indices) =
        generate_from_pixel_data(data, width, height, threshold, density, border_px)?;
    params
        .push_result_array_float(&vertices)
        .context("failed to return mesh vertices")?;
    params
        .push_result_array_int(&indices)
        .context("failed to return mesh indices")?;
    Ok(())
}

fn set_callback_error(raw: *mut SCRIPT_MODULE_PARAM, message: &str) {
    if raw.is_null() {
        return;
    }

    // CString cannot contain embedded NULs. Replacing them still gives the Lua
    // caller a useful error instead of panicking while reporting an error.
    let message = message.replace('\0', "\\0");
    if let Ok(message) = CString::new(message) {
        // SAFETY: This is only called synchronously from the host callback, so
        // the SDK table and the temporary C string are both valid for the call.
        unsafe { ((*raw).set_error)(message.as_ptr()) };
    }
}

fn generate_from_pixel_data(
    data: NonNull<u8>,
    width: i32,
    height: i32,
    threshold: i32,
    density: i32,
    border_px: i32,
) -> AnyResult<(Vec<f64>, Vec<i32>)> {
    let (width, height, byte_len) = checked_image_layout(width, height)?;

    // SAFETY: AviUtl2 owns this buffer. `obj.getpixeldata()` returns a pointer
    // to width * height RGBA32 pixels, and Lua passes that pointer and the same
    // dimensions in this synchronous call. The slice is never retained.
    let rgba = unsafe { std::slice::from_raw_parts(data.as_ptr(), byte_len) };
    let mesh = generate_mesh(rgba, width, height, threshold, density, border_px)?;

    let mut vertices = Vec::with_capacity(mesh.points.len() * 2);
    for (x, y) in mesh.points {
        vertices.push(f64::from(x));
        vertices.push(f64::from(y));
    }

    let mut indices = Vec::with_capacity(mesh.triangles.len() * 3);
    for triangle in mesh.triangles {
        for index in triangle {
            indices.push(i32::try_from(index).context("mesh has too many vertices for AviUtl2")?);
        }
    }

    Ok((vertices, indices))
}

struct Mesh {
    points: Vec<(f32, f32)>,
    triangles: Vec<[usize; 3]>,
}

fn checked_image_layout(width: i32, height: i32) -> AnyResult<(usize, usize, usize)> {
    anyhow::ensure!(width > 0 && height > 0, "image dimensions must be positive");

    let width = usize::try_from(width).context("invalid image width")?;
    let height = usize::try_from(height).context("invalid image height")?;
    let byte_len = width
        .checked_mul(height)
        .and_then(|pixels| pixels.checked_mul(RGBA_CHANNELS))
        .context("RGBA image dimensions are too large")?;

    // Rust slices use isize-sized offsets. Reject impossible layouts before
    // constructing one from the host-owned pointer.
    anyhow::ensure!(byte_len <= isize::MAX as usize, "RGBA image is too large");
    Ok((width, height, byte_len))
}

fn generate_mesh(
    rgba: &[u8],
    width: usize,
    height: usize,
    threshold: i32,
    density: i32,
    border_px: i32,
) -> AnyResult<Mesh> {
    let expected_len = width
        .checked_mul(height)
        .and_then(|pixels| pixels.checked_mul(RGBA_CHANNELS))
        .context("RGBA image dimensions are too large")?;
    anyhow::ensure!(
        rgba.len() == expected_len,
        "RGBA buffer length does not match its dimensions"
    );

    let threshold = threshold.clamp(0, 255) as u8;
    let density = density.clamp(MIN_DENSITY, MAX_DENSITY) as usize;
    let border = border_px.max(0) as f32;

    let max_dim = width.max(height) as f32;
    let min_spacing = max_dim / density as f32;
    let epsilon = min_spacing * 0.4;

    // Douglas-Peucker may move a contour inward by epsilon. Dilating by the
    // requested border plus epsilon prevents that simplification from cutting
    // into opaque pixels.
    let total_dilate = border + epsilon;
    let pad = total_dilate.ceil() as usize + 2;
    let padded_width = width
        .checked_add(pad.checked_mul(2).context("mesh padding is too large")?)
        .context("padded image width is too large")?;
    let padded_height = height
        .checked_add(pad.checked_mul(2).context("mesh padding is too large")?)
        .context("padded image height is too large")?;
    let padded_len = padded_width
        .checked_mul(padded_height)
        .context("padded image is too large")?;

    let mut binary = vec![false; padded_len];
    for y in 0..height {
        let src_row = y * width * RGBA_CHANNELS;
        let dst_row = (y + pad) * padded_width + pad;
        for x in 0..width {
            binary[dst_row + x] = rgba[src_row + x * RGBA_CHANNELS + 3] >= threshold;
        }
    }

    if total_dilate > 0.0 {
        dilate::dilate_edt(&mut binary, padded_width, padded_height, total_dilate);
    }

    let dilated_alpha: Vec<u8> = binary
        .iter()
        .map(|&opaque| if opaque { 255 } else { 0 })
        .collect();
    let contours = contour::extract_contours_from_binary(&binary, padded_width, padded_height);

    let simplified_contours: Vec<Vec<(f32, f32)>> = contours
        .iter()
        .filter_map(|points| {
            let simplified = if points.len() > 3 {
                simplify::simplify_loop(points, epsilon)
            } else {
                points.clone()
            };
            (simplified.len() >= 3).then_some(simplified)
        })
        .collect();

    let mut points = Vec::new();
    let mut constraint_edges = Vec::new();
    for contour in simplified_contours {
        let base = points.len();
        let count = contour.len();
        points.extend(contour);
        for index in 0..count {
            constraint_edges.push((base + index, base + (index + 1) % count));
        }
    }

    let seed = (width as u64 * 31_337 + height as u64 * 7_919 + density as u64 * 104_729) | 1;
    points.extend(poisson::poisson_disk_sample(
        padded_width,
        padded_height,
        min_spacing,
        &dilated_alpha,
        1,
        &points,
        seed,
    ));

    if points.len() < 3 {
        return Ok(Mesh {
            points: Vec::new(),
            triangles: Vec::new(),
        });
    }

    let triangles = delaunay::triangulate(
        &points,
        &constraint_edges,
        &dilated_alpha,
        padded_width,
        padded_height,
        1,
    );

    if triangles.is_empty() {
        points.clear();
    }

    Ok(Mesh { points, triangles })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_pipeline_generates_valid_indices() {
        let width = 32;
        let height = 32;
        let mut rgba = vec![0_u8; width * height * RGBA_CHANNELS];
        for y in 0..height {
            for x in 0..width {
                let dx = x as f32 - width as f32 * 0.5;
                let dy = y as f32 - height as f32 * 0.5;
                if dx * dx + dy * dy <= 12.0 * 12.0 {
                    rgba[(y * width + x) * RGBA_CHANNELS + 3] = 255;
                }
            }
        }

        let mesh = generate_mesh(&rgba, width, height, 128, 10, 0).unwrap();
        assert!(!mesh.points.is_empty());
        assert!(!mesh.triangles.is_empty());
        assert!(
            mesh.triangles
                .iter()
                .flatten()
                .all(|&index| index < mesh.points.len())
        );
    }

    #[test]
    fn transparent_image_returns_an_empty_mesh() {
        let rgba = vec![0_u8; 16 * 16 * RGBA_CHANNELS];
        let mesh = generate_mesh(&rgba, 16, 16, 128, 10, 0).unwrap();
        assert!(mesh.points.is_empty());
        assert!(mesh.triangles.is_empty());
    }

    #[test]
    fn rejects_mismatched_buffer_length() {
        let error = generate_mesh(&[0; 3], 1, 1, 128, 10, 0)
            .err()
            .expect("invalid buffer should fail");
        assert!(error.to_string().contains("buffer length"));
    }

    #[test]
    fn rejects_invalid_dimensions() {
        assert!(checked_image_layout(0, 10).is_err());
        assert!(checked_image_layout(10, -1).is_err());
    }
}
