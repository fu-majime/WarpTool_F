//! Shared user frame input for ARAP constraints and rigid MLS observations.
//! Position and frame constraints remain separate, leaving bend translations free.
use aviutl2::anyhow;
pub const STRIDE: usize = 12;
#[derive(Clone, Copy)]
pub struct Pose {
    pub source: (f64, f64),
    pub center: (f64, f64),
    pub natural_angle: f64,
    pub angle: f64,
    pub scale: f64,
    pub radius: f64,
    pub moving: bool,
    pub rotate: bool,
    pub resize: bool,
}
pub fn parse(values: &[f64]) -> anyhow::Result<Vec<Pose>> {
    anyhow::ensure!(
        values.len().is_multiple_of(STRIDE) && values.iter().all(|v| v.is_finite()),
        "invalid pose payload"
    );
    values
        .chunks_exact(STRIDE)
        .map(|v| {
            anyhow::ensure!(
                v[5] >= 0.0 && v[7] >= 0.0 && v[8] > 0.0,
                "invalid pose scale or radius"
            );
            Ok(Pose {
                source: (v[0], v[1]),
                center: (v[2], v[3]),
                natural_angle: v[4],
                angle: v[6],
                scale: v[7],
                radius: v[8],
                moving: v[9] != 0.0,
                rotate: v[10] != 0.0,
                resize: v[11] != 0.0,
            })
        })
        .collect()
}
pub fn weight(distance: f64, radius: f64) -> f64 {
    let t = ((2.0 * radius - distance) / radius).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}
