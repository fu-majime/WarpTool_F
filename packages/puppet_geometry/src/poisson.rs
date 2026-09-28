//! Deterministic variable-radius sampling, dense near boundaries and sparse
//! in broad interiors. Scan candidates cover every disconnected component.
pub fn poisson_disk_sample(
    w: usize,
    h: usize,
    base: f32,
    alpha: &[u8],
    threshold: u8,
    contour: &[(f32, f32)],
    seed: u64,
) -> Vec<(f32, f32)> {
    let base = base.max(1.0);
    let mut depth = vec![0.0_f32; w * h];
    for i in 0..w * h {
        if alpha[i] >= threshold {
            depth[i] = (w + h) as f32;
        }
    }
    for y in 0..h {
        for x in 0..w {
            let i = y * w + x;
            if x > 0 {
                depth[i] = depth[i].min(depth[i - 1] + 1.0);
            }
            if y > 0 {
                depth[i] = depth[i].min(depth[i - w] + 1.0);
            }
        }
    }
    for y in (0..h).rev() {
        for x in (0..w).rev() {
            let i = y * w + x;
            if x + 1 < w {
                depth[i] = depth[i].min(depth[i + 1] + 1.0);
            }
            if y + 1 < h {
                depth[i] = depth[i].min(depth[i + w] + 1.0);
            }
        }
    }
    let cell = base * 0.5;
    let gw = (w as f32 / cell).ceil() as usize + 1;
    let gh = (h as f32 / cell).ceil() as usize + 1;
    let mut grid = vec![Vec::<(f32, f32, f32)>::new(); gw * gh];
    let (hw, hh) = (w as f32 * 0.5, h as f32 * 0.5);
    for &(x, y) in contour {
        let (x, y) = (x + hw, y + hh);
        let (gx, gy) = ((x / cell) as usize, (y / cell) as usize);
        if gx < gw && gy < gh {
            grid[gy * gw + gx].push((x, y, 0.0));
        }
    }
    let step = (base * 0.4).floor().max(1.0) as usize;
    let mut candidates = Vec::new();
    let mut random = seed.max(1);
    let mut rand = || {
        random ^= random << 13;
        random ^= random >> 7;
        random ^= random << 17;
        random as f64 / u64::MAX as f64
    };
    for y in (0..h).step_by(step) {
        for x in (0..w).step_by(step) {
            let px = (x as f64 + 0.5 + rand() * (step - 1) as f64) as f32;
            let py = (y as f64 + 0.5 + rand() * (step - 1) as f64) as f32;
            if px >= w as f32 || py >= h as f32 {
                continue;
            }
            let i = py as usize * w + px as usize;
            if alpha[i] < threshold {
                continue;
            }
            let radius = base * (1.0 + (depth[i] / base).min(1.0));
            candidates.push((px, py, radius));
        }
    }
    candidates.sort_by(|a, b| a.2.total_cmp(&b.2));
    let mut result = Vec::new();
    for (x, y, radius) in candidates {
        let (gx, gy) = ((x / cell) as usize, (y / cell) as usize);
        let mut valid = true;
        'neighbors: for yy in gy.saturating_sub(4)..=(gy + 4).min(gh - 1) {
            for xx in gx.saturating_sub(4)..=(gx + 4).min(gw - 1) {
                for &(px, py, pr) in &grid[yy * gw + xx] {
                    let spacing = if pr == 0.0 {
                        radius * 0.5
                    } else {
                        radius.max(pr)
                    };
                    if (x - px).powi(2) + (y - py).powi(2) < spacing * spacing {
                        valid = false;
                        break 'neighbors;
                    }
                }
            }
        }
        if valid {
            grid[gy * gw + gx].push((x, y, radius));
            result.push((x - hw, y - hh));
        }
    }
    result
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn samples_all_islands_and_stays_inside() {
        let (w, h) = (160, 80);
        let mut alpha = vec![0; w * h];
        for y in 5..75 {
            for x in 5..65 {
                alpha[y * w + x] = 255;
                alpha[y * w + x + 90] = 255;
            }
        }
        let p = poisson_disk_sample(w, h, 8.0, &alpha, 1, &[], 42);
        assert!(p.iter().any(|p| p.0 < 0.0) && p.iter().any(|p| p.0 > 0.0));
        for &(x, y) in &p {
            assert_eq!(alpha[(y + 40.0) as usize * w + (x + 80.0) as usize], 255);
        }
        assert!(p.len() < 160);
        assert_eq!(p, poisson_disk_sample(w, h, 8.0, &alpha, 1, &[], 42));
    }
    #[test]
    fn empty_has_no_samples() {
        assert!(poisson_disk_sample(16, 16, 2.0, &[0; 256], 1, &[], 1).is_empty());
    }
}
