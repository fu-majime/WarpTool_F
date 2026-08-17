// 入力画像キャプチャ vramwiz氏を参考
use std::sync::Mutex;

use aviutl2::filter::FilterProcVideo;
use half::f16;
use windows::{
    Win32::Graphics::{
        Direct3D11::{
            D3D11_CPU_ACCESS_READ, D3D11_MAP_READ, D3D11_MAPPED_SUBRESOURCE, D3D11_TEXTURE2D_DESC,
            D3D11_USAGE_STAGING, ID3D11Device, ID3D11Texture2D,
        },
        Dxgi::Common::{
            DXGI_FORMAT, DXGI_FORMAT_B8G8R8A8_UNORM, DXGI_FORMAT_B8G8R8A8_UNORM_SRGB,
            DXGI_FORMAT_R8G8B8A8_UNORM, DXGI_FORMAT_R8G8B8A8_UNORM_SRGB,
            DXGI_FORMAT_R16G16B16A16_FLOAT, DXGI_FORMAT_R16G16B16A16_UNORM,
        },
    },
    core::Interface,
};

use crate::state::{CaptureSource, SharedState};

const MAX_CAPTURE_DIMENSION: u32 = 16_384;

pub struct FrameCapture {
    resources: Mutex<CaptureResources>,
}

struct CaptureResources {
    staging: Option<ID3D11Texture2D>,
    device_identity: usize,
    width: u32,
    height: u32,
    format: DXGI_FORMAT,
}

impl FrameCapture {
    pub fn new() -> Self {
        Self {
            resources: Mutex::new(CaptureResources {
                staging: None,
                device_identity: 0,
                width: 0,
                height: 0,
                format: DXGI_FORMAT(0),
            }),
        }
    }

    pub fn capture<U: aviutl2::filter::FilterUserdata>(
        &self,
        video: &mut FilterProcVideo<U>,
        shared: &SharedState,
    ) {
        let input = video.get_image_texture2d();
        let (raw, source) = if input.is_null() {
            (
                video.get_framebuffer_texture2d(),
                CaptureSource::FramebufferFallback,
            )
        } else {
            (input, CaptureSource::InputImage)
        };
        if raw.is_null() {
            shared.publish_capture_error("入力画像テクスチャを取得できませんでした。");
            return;
        }

        let result = unsafe {
            // AviUtl2 owns this COM reference and guarantees it only for proc_video.
            // A borrowed wrapper prevents Release from being called on that reference.
            let Some(texture) = ID3D11Texture2D::from_raw_borrowed(&raw) else {
                return shared.publish_capture_error("入力画像テクスチャが無効です。");
            };
            self.copy_texture(texture)
        };
        match result {
            Ok((width, height, pixels)) => shared.publish_capture(width, height, pixels, source),
            Err(error) => shared.publish_capture_error(error),
        }
    }

    unsafe fn copy_texture(
        &self,
        source: &ID3D11Texture2D,
    ) -> Result<(usize, usize, Vec<u8>), String> {
        let mut desc = D3D11_TEXTURE2D_DESC::default();
        unsafe { source.GetDesc(&mut desc) };
        validate_description(&desc)?;
        let pixel_format = PixelFormat::try_from(desc.Format)?;
        let row_bytes = desc.Width as usize * pixel_format.bytes_per_pixel();

        let device = unsafe { source.GetDevice() }
            .map_err(|error| format!("D3D11デバイスを取得できませんでした: {error}"))?;
        let device_identity = Interface::as_raw(&device) as usize;

        let mut resources = self.resources.lock().unwrap();
        let needs_recreate = resources.staging.is_none()
            || resources.device_identity != device_identity
            || resources.width != desc.Width
            || resources.height != desc.Height
            || resources.format != desc.Format;
        if needs_recreate {
            resources.staging = Some(unsafe { create_staging_texture(&device, &desc)? });
            resources.device_identity = device_identity;
            resources.width = desc.Width;
            resources.height = desc.Height;
            resources.format = desc.Format;
        }
        let staging = resources.staging.as_ref().unwrap();
        let context = unsafe { device.GetImmediateContext() }
            .map_err(|error| format!("D3D11コンテキストを取得できませんでした: {error}"))?;
        unsafe { context.CopyResource(staging, source) };

        let mut mapped = D3D11_MAPPED_SUBRESOURCE::default();
        unsafe { context.Map(staging, 0, D3D11_MAP_READ, 0, Some(&mut mapped)) }
            .map_err(|error| format!("入力画像をMapできませんでした: {error}"))?;

        let result = (|| {
            let row_pitch = mapped.RowPitch as usize;
            if mapped.pData.is_null() || row_pitch < row_bytes {
                return Err("入力画像のRowPitchが不正です。".to_owned());
            }
            let data_len = row_pitch
                .checked_mul(desc.Height as usize)
                .ok_or_else(|| "入力画像のサイズが大きすぎます。".to_owned())?;
            let data = unsafe { std::slice::from_raw_parts(mapped.pData.cast::<u8>(), data_len) };
            convert_rows(
                data,
                row_pitch,
                desc.Width as usize,
                desc.Height as usize,
                pixel_format,
            )
        })();
        unsafe { context.Unmap(staging, 0) };
        result.map(|pixels| (desc.Width as usize, desc.Height as usize, pixels))
    }
}

fn validate_description(desc: &D3D11_TEXTURE2D_DESC) -> Result<(), String> {
    if desc.Width == 0
        || desc.Height == 0
        || desc.Width > MAX_CAPTURE_DIMENSION
        || desc.Height > MAX_CAPTURE_DIMENSION
    {
        return Err(format!(
            "入力画像の寸法が不正です: {} x {}",
            desc.Width, desc.Height
        ));
    }
    if desc.SampleDesc.Count != 1 {
        return Err(format!(
            "MSAA入力画像には対応していません: {} samples",
            desc.SampleDesc.Count
        ));
    }
    Ok(())
}

unsafe fn create_staging_texture(
    device: &ID3D11Device,
    source: &D3D11_TEXTURE2D_DESC,
) -> Result<ID3D11Texture2D, String> {
    let mut desc = *source;
    desc.MipLevels = 1;
    desc.ArraySize = 1;
    desc.Usage = D3D11_USAGE_STAGING;
    desc.BindFlags = 0;
    desc.CPUAccessFlags = D3D11_CPU_ACCESS_READ.0 as u32;
    desc.MiscFlags = 0;
    let mut staging = None;
    unsafe { device.CreateTexture2D(&desc, None, Some(&mut staging)) }
        .map_err(|error| format!("staging textureを作成できませんでした: {error}"))?;
    staging.ok_or_else(|| "staging textureが返されませんでした。".to_owned())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PixelFormat {
    Rgba8,
    Bgra8,
    Rgba16Unorm,
    Rgba16Float,
}

impl PixelFormat {
    fn bytes_per_pixel(self) -> usize {
        match self {
            Self::Rgba8 | Self::Bgra8 => 4,
            Self::Rgba16Unorm | Self::Rgba16Float => 8,
        }
    }
}

impl TryFrom<DXGI_FORMAT> for PixelFormat {
    type Error = String;

    fn try_from(value: DXGI_FORMAT) -> Result<Self, Self::Error> {
        match value {
            DXGI_FORMAT_R8G8B8A8_UNORM | DXGI_FORMAT_R8G8B8A8_UNORM_SRGB => Ok(Self::Rgba8),
            DXGI_FORMAT_B8G8R8A8_UNORM | DXGI_FORMAT_B8G8R8A8_UNORM_SRGB => Ok(Self::Bgra8),
            DXGI_FORMAT_R16G16B16A16_UNORM => Ok(Self::Rgba16Unorm),
            DXGI_FORMAT_R16G16B16A16_FLOAT => Ok(Self::Rgba16Float),
            _ => Err(format!("未対応のDXGI形式です: {}", value.0)),
        }
    }
}

fn convert_rows(
    source: &[u8],
    row_pitch: usize,
    width: usize,
    height: usize,
    format: PixelFormat,
) -> Result<Vec<u8>, String> {
    let source_row_bytes = width
        .checked_mul(format.bytes_per_pixel())
        .ok_or_else(|| "入力画像の行サイズが大きすぎます。".to_owned())?;
    if row_pitch < source_row_bytes || source.len() < row_pitch.saturating_mul(height) {
        return Err("入力画像バッファが不足しています。".to_owned());
    }
    let output_len = width
        .checked_mul(height)
        .and_then(|value| value.checked_mul(4))
        .ok_or_else(|| "入力画像が大きすぎます。".to_owned())?;
    let mut output = Vec::with_capacity(output_len);

    for y in 0..height {
        let row = &source[y * row_pitch..y * row_pitch + source_row_bytes];
        match format {
            PixelFormat::Rgba8 => output.extend_from_slice(row),
            PixelFormat::Bgra8 => {
                for pixel in row.chunks_exact(4) {
                    output.extend_from_slice(&[pixel[2], pixel[1], pixel[0], pixel[3]]);
                }
            }
            PixelFormat::Rgba16Unorm => {
                for pixel in row.chunks_exact(8) {
                    for channel in pixel.chunks_exact(2) {
                        let value = u16::from_le_bytes([channel[0], channel[1]]);
                        output.push(((u32::from(value) * 255 + 32_767) / 65_535) as u8);
                    }
                }
            }
            PixelFormat::Rgba16Float => {
                for pixel in row.chunks_exact(8) {
                    for channel in pixel.chunks_exact(2) {
                        let value = f16::from_bits(u16::from_le_bytes([channel[0], channel[1]]));
                        output.push(float_to_byte(value.to_f32()));
                    }
                }
            }
        }
    }
    Ok(output)
}

fn float_to_byte(value: f32) -> u8 {
    if !value.is_finite() || value <= 0.0 {
        0
    } else if value >= 1.0 {
        255
    } else {
        (value * 255.0).round() as u8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn converts_bgra_with_row_padding() {
        let source = [1, 2, 3, 4, 0xaa, 0xbb, 0xcc, 0xdd, 5, 6, 7, 8, 0, 0, 0, 0];
        let output = convert_rows(&source, 8, 1, 2, PixelFormat::Bgra8).unwrap();
        assert_eq!(output, [3, 2, 1, 4, 7, 6, 5, 8]);
    }

    #[test]
    fn copies_rgba_while_discarding_row_padding() {
        let source = [1, 2, 3, 4, 9, 9, 9, 9, 5, 6, 7, 8, 9, 9, 9, 9];
        let output = convert_rows(&source, 8, 1, 2, PixelFormat::Rgba8).unwrap();
        assert_eq!(output, [1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn converts_unorm_and_half_float() {
        let unorm = [0, 0, 0xff, 0xff, 0x80, 0x80, 0, 0];
        assert_eq!(
            convert_rows(&unorm, 8, 1, 1, PixelFormat::Rgba16Unorm).unwrap(),
            [0, 255, 128, 0]
        );

        let mut half = Vec::new();
        for value in [-1.0, 0.5, 1.0, f32::NAN] {
            half.extend_from_slice(&f16::from_f32(value).to_bits().to_le_bytes());
        }
        assert_eq!(
            convert_rows(&half, 8, 1, 1, PixelFormat::Rgba16Float).unwrap(),
            [0, 128, 255, 0]
        );
    }
}
