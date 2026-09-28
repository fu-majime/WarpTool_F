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

mod arap;
mod contour;
mod deformation;
mod delaunay;
mod dilate;
mod embedding;
mod linear;
mod math;
mod mls;
mod poisson;
mod pose;
#[cfg(test)]
mod regression;
mod simplify;

const RGBA_CHANNELS: usize = 4;
const MIN_DENSITY: i32 = 5;
const MAX_DENSITY: i32 = 200;

#[aviutl2::plugin(ScriptModule)]
struct PuppetGeometryModule;

impl aviutl2::module::ScriptModule for PuppetGeometryModule {
    fn new(_info: aviutl2::AviUtl2Info) -> AnyResult<Self> {
        Ok(Self)
    }

    fn plugin_info(&self) -> aviutl2::module::ScriptModuleTable {
        aviutl2::module::ScriptModuleTable {
            information: format!(
                "Puppet geometry engine for AviUtl2 / v{} / Fu-Majime",
                env!("CARGO_PKG_VERSION")
            ),
            functions: Self::functions(),
        }
    }
}

impl ScriptModuleFunctions for PuppetGeometryModule {
    fn functions() -> Vec<ModuleFunction> {
        vec![
            ModuleFunction {
                name: "prepare_mesh".to_owned(),
                func: prepare_callback,
            },
            ModuleFunction {
                name: "generate".to_owned(),
                func: generate_callback,
            },
            ModuleFunction {
                name: "deform_mls".to_owned(),
                func: deform_callback,
            },
            ModuleFunction {
                name: "deform_arap".to_owned(),
                func: deform_arap_callback,
            },
        ]
    }
}

extern "C" fn prepare_callback(raw: *mut SCRIPT_MODULE_PARAM) {
    let result = catch_unwind(AssertUnwindSafe(|| -> AnyResult<()> {
        anyhow::ensure!(!raw.is_null(), "null module-call handle");
        // SAFETY: synchronous host-owned call table, as in deform_callback.
        let mut params = unsafe { ScriptModuleCallHandle::from_raw(raw) };
        let floats = |index| {
            (0..params.get_param_array_len(index))
                .map(|i| params.get_param_array_float(index, i))
                .collect::<Vec<_>>()
        };
        let vertices = floats(0);
        let pins = floats(2);
        let indices = (0..params.get_param_array_len(1))
            .map(|i| params.get_param_array_int(1, i))
            .collect::<Vec<_>>();
        let radius = params
            .get_param_float(3)
            .context("invalid refinement radius")?;
        let spacing = params.get_param_float(4).unwrap_or(0.0);
        let (vertices, indices) =
            embedding::prepare_with_spacing(&vertices, &indices, &pins, radius, spacing)?;
        params.push_result_array_float(&vertices)?;
        params.push_result_array_int(&indices)?;
        Ok(())
    }));
    match result {
        Ok(Ok(())) => {}
        Ok(Err(e)) => set_callback_error(raw, &e.to_string()),
        Err(_) => set_callback_error(raw, "internal panic in prepare_mesh"),
    }
}

extern "C" fn deform_arap_callback(raw: *mut SCRIPT_MODULE_PARAM) {
    let result = catch_unwind(AssertUnwindSafe(|| deform_callback_inner(raw, true)));
    match result {
        Ok(Ok(())) => {}
        Ok(Err(error)) => set_callback_error(raw, &error.to_string()),
        Err(_) => set_callback_error(raw, "internal panic in deform_arap"),
    }
}

extern "C" fn deform_callback(raw: *mut SCRIPT_MODULE_PARAM) {
    let result = catch_unwind(AssertUnwindSafe(|| deform_callback_inner(raw, false)));
    match result {
        Ok(Ok(())) => {}
        Ok(Err(error)) => set_callback_error(raw, &error.to_string()),
        Err(_) => set_callback_error(raw, "internal panic in deform_mls"),
    }
}

fn deform_callback_inner(raw: *mut SCRIPT_MODULE_PARAM, arap: bool) -> AnyResult<()> {
    anyhow::ensure!(!raw.is_null(), "AviUtl2 passed a null module-call handle");
    // SAFETY: the host owns the call table for the duration of this callback.
    let mut params = unsafe { ScriptModuleCallHandle::from_raw(raw) };
    let floats = |index| {
        (0..params.get_param_array_len(index))
            .map(|i| params.get_param_array_float(index, i))
            .collect::<Vec<_>>()
    };
    let ints = |index| {
        (0..params.get_param_array_len(index))
            .map(|i| params.get_param_array_int(index, i))
            .collect::<Vec<_>>()
    };
    let vertices = floats(0);
    let indices = ints(1);
    let sources = floats(2);
    let destinations = floats(3);
    let layers = floats(4);
    let poses = floats(9);
    let divisions = params
        .get_param_int(6)
        .context("invalid subdivision count")?;
    let width = params.get_param_float(7).context("invalid width")?;
    let height = params.get_param_float(8).context("invalid height")?;
    let method = if arap {
        deformation::Method::Arap
    } else {
        deformation::Method::RigidMls(params.get_param_float(5).context("invalid stiffness")?)
    };
    let deformation = deformation::deform(
        &vertices,
        &indices,
        &sources,
        &destinations,
        &layers,
        divisions,
        width,
        height,
        &poses,
        method,
    )?;
    params
        .push_result_array_float(&deformation.vertices)
        .context("failed to return deformed vertices")?;
    params
        .push_result_array_float(&deformation.render_vertices)
        .context("failed to return render vertices")?;
    params
        .push_result_array_float(&deformation.wire_vertices)
        .context("failed to return wire vertices")?;
    Ok(())
}

aviutl2::register_script_module!(PuppetGeometryModule);

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

#[derive(Clone)]
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
    #[derive(PartialEq, Eq)]
    struct Key {
        width: usize,
        height: usize,
        density: i32,
        border: i32,
        mask: Vec<u64>,
    }
    // Equal binary masks generate equal geometry, even at different cutoffs.
    thread_local! { static CACHE:std::cell::RefCell<std::collections::VecDeque<(Key,Mesh)>>=const { std::cell::RefCell::new(std::collections::VecDeque::new()) }; }
    let pixels = width.checked_mul(height).context("image too large")?;
    anyhow::ensure!(
        pixels.checked_mul(4) == Some(rgba.len()),
        "RGBA buffer length does not match its dimensions"
    );
    let mut mask = vec![0_u64; pixels.div_ceil(64)];
    let cutoff = threshold.clamp(0, 255) as u8;
    for (i, pixel) in rgba.chunks_exact(4).enumerate() {
        if pixel[3] >= cutoff {
            mask[i / 64] |= 1_u64 << (i % 64);
        }
    }
    let key = Key {
        width,
        height,
        density: density.clamp(MIN_DENSITY, MAX_DENSITY),
        border: border_px.max(0),
        mask,
    };
    if let Some(mesh) = CACHE.with(|cache| {
        let mut cache = cache.borrow_mut();
        let i = cache.iter().position(|(k, _)| k == &key)?;
        let entry = cache.remove(i)?;
        let result = entry.1.clone();
        cache.push_front(entry);
        Some(result)
    }) {
        return Ok(mesh);
    }
    let mesh = generate_mesh_uncached(rgba, width, height, threshold, density, border_px)?;
    CACHE.with(|cache| {
        let mut cache = cache.borrow_mut();
        cache.push_front((key, mesh.clone()));
        cache.truncate(4);
    });
    Ok(mesh)
}

fn generate_mesh_uncached(
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
    let epsilon = (min_spacing * 0.2).max(1.25);
    let area_budget = 0.25 / density as f64;
    // Small enclosed holes remain visible in the alpha texture without forcing
    // sub-spacing rings into a coarse deformation mesh. Open notches stay open.
    let hole_span = (min_spacing * 0.125 * ((10.0 - density as f32) / 5.0).clamp(0.0, 1.0)).floor();

    // Only the explicit border expands the outer silhouette.
    let total_dilate = border;
    let pad = (total_dilate + hole_span).ceil() as usize + 2;
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
            let alpha = rgba[src_row + x * RGBA_CHANNELS + 3];
            binary[dst_row + x] = alpha >= threshold;
        }
    }

    if total_dilate > 0.0 {
        dilate::dilate_edt(&mut binary, padded_width, padded_height, total_dilate);
    }
    if hole_span >= 1.0 {
        let mut visited = vec![false; binary.len()];
        for seed in 0..binary.len() {
            if binary[seed] || visited[seed] {
                continue;
            }
            let mut cells = vec![seed];
            visited[seed] = true;
            let mut cursor = 0;
            let mut outside = false;
            while cursor < cells.len() {
                let i = cells[cursor];
                cursor += 1;
                let (x, y) = (i % padded_width, i / padded_width);
                if x == 0 || y == 0 || x + 1 == padded_width || y + 1 == padded_height {
                    outside = true;
                }
                for j in [
                    i.checked_sub(1).filter(|_| x > 0),
                    i.checked_add(1).filter(|_| x + 1 < padded_width),
                    i.checked_sub(padded_width),
                    i.checked_add(padded_width).filter(|&j| j < binary.len()),
                ]
                .into_iter()
                .flatten()
                {
                    if !binary[j] && !visited[j] {
                        visited[j] = true;
                        cells.push(j);
                    }
                }
            }
            if !outside && cells.len() as f32 <= hole_span * hole_span {
                for i in cells {
                    binary[i] = true;
                }
            }
        }
    }

    let dilated_alpha: Vec<u8> = binary
        .iter()
        .map(|&opaque| if opaque { 255 } else { 0 })
        .collect();
    let contours = contour::extract_contours_from_binary(&binary, padded_width, padded_height);

    let simplified_contours = delaunay::validated_contours(
        &contours,
        simplify::contours(&contours, epsilon, min_spacing * 1.5, area_budget),
    )?;
    let seed = (width as u64 * 31_337 + height as u64 * 7_919 + density as u64 * 104_729) | 1;
    let (mut points, triangles) = {
        let mut points = Vec::new();
        let mut constraint_edges = Vec::new();
        for contour in &simplified_contours {
            let base = points.len();
            let count = contour.len();
            points.extend(contour);
            for i in 0..count {
                constraint_edges.push((base + i, base + (i + 1) % count));
            }
        }
        let interior = poisson::poisson_disk_sample(
            padded_width,
            padded_height,
            min_spacing,
            &dilated_alpha,
            1,
            &points,
            seed,
        );
        // Boundary vertices alone do not provide clearance from a long edge.
        // Reject near-edge samples that would create almost-flat sliver faces.
        let interior = interior
            .into_iter()
            .filter(|&(x, y)| {
                constraint_edges.iter().all(|&(a, b)| {
                    let a = points[a];
                    let b = points[b];
                    let (u, v) = (b.0 - a.0, b.1 - a.1);
                    let t = (((x - a.0) * u + (y - a.1) * v) / (u * u + v * v).max(1e-12))
                        .clamp(0.0, 1.0);
                    (x - a.0 - t * u).hypot(y - a.1 - t * v) >= min_spacing * 0.25
                })
            })
            .collect::<Vec<_>>();
        points.extend(interior);
        if points.len() < 3 {
            return Ok(Mesh {
                points: Vec::new(),
                triangles: Vec::new(),
            });
        }
        let triangles = delaunay::triangulate(&points, &constraint_edges)?;
        (points, triangles)
    };

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
