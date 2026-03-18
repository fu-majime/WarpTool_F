/// Poisson Disk Sampling (Bridsonのアルゴリズム)
///
/// 均一で自然なランダム点分布を生成する。
/// 格子状にならないのが格子サンプリングとの最大の違い。

/// Bridsonの高速Poisson Disk Sampling。
///
/// - `w`, `h`: サンプリング領域サイズ（ピクセル）
/// - `min_dist`: 点間の最小距離
/// - `alpha`: α画像データ (w*h), しきい値以上の領域にのみサンプリング
/// - `threshold`: αしきい値
/// - `contour_pts`: 輪郭上の点（これらとの距離も min_dist を保つ）
/// - `seed`: 乱数シード
///
/// 戻り値: サンプリングされた内部点のリスト（画像中心原点）
pub fn poisson_disk_sample(
    w: usize,
    h: usize,
    min_dist: f32,
    alpha: &[u8],
    threshold: u8,
    contour_pts: &[(f32, f32)],
    seed: u64,
) -> Vec<(f32, f32)> {
    let hw = w as f32 * 0.5;
    let hh = h as f32 * 0.5;

    // 加速グリッド
    let cell_size = min_dist / std::f32::consts::SQRT_2;
    let grid_w = (w as f32 / cell_size).ceil() as usize + 1;
    let grid_h = (h as f32 / cell_size).ceil() as usize + 1;
    let mut grid: Vec<Option<usize>> = vec![None; grid_w * grid_h];

    // 結果リスト
    let mut points: Vec<(f32, f32)> = Vec::new();
    // アクティブリスト
    let mut active: Vec<usize> = Vec::new();

    let max_attempts = 30;
    let min_dist2 = min_dist * min_dist;

    // 簡易PRNG (xorshift64)
    let mut rng_state = if seed == 0 { 0x12345678u64 } else { seed };
    let mut next_rand = move || -> f64 {
        rng_state ^= rng_state << 13;
        rng_state ^= rng_state >> 7;
        rng_state ^= rng_state << 17;
        (rng_state as f64) / (u64::MAX as f64)
    };

    // 輪郭点をグリッドに登録
    let mut all_pts: Vec<(f32, f32)> = Vec::new();
    for &(cx, cy) in contour_pts {
        // 画像座標に変換 (中心原点 → 左上原点)
        let px = cx + hw;
        let py = cy + hh;
        if px >= 0.0 && px < w as f32 && py >= 0.0 && py < h as f32 {
            let idx = all_pts.len();
            all_pts.push((px, py));
            let gx = (px / cell_size) as usize;
            let gy = (py / cell_size) as usize;
            if gx < grid_w && gy < grid_h {
                grid[gy * grid_w + gx] = Some(idx);
            }
        }
    }

    // α領域内のランダムな初期点
    let mut initial_found = false;
    for _ in 0..1000 {
        let px = next_rand() as f32 * w as f32;
        let py = next_rand() as f32 * h as f32;
        let ix = px as usize;
        let iy = py as usize;
        if ix < w && iy < h && alpha[iy * w + ix] >= threshold {
            // 輪郭点との距離チェック
            if is_valid_point(px, py, &all_pts, &grid, grid_w, grid_h, cell_size, min_dist2) {
                let idx = all_pts.len();
                all_pts.push((px, py));
                let gx = (px / cell_size) as usize;
                let gy = (py / cell_size) as usize;
                if gx < grid_w && gy < grid_h {
                    grid[gy * grid_w + gx] = Some(idx);
                }
                points.push((px, py));
                active.push(idx);
                initial_found = true;
                break;
            }
        }
    }

    if !initial_found {
        // α領域の中心付近を強制使用
        let cx = w as f32 * 0.5;
        let cy = h as f32 * 0.5;
        let idx = all_pts.len();
        all_pts.push((cx, cy));
        let gx = (cx / cell_size) as usize;
        let gy = (cy / cell_size) as usize;
        if gx < grid_w && gy < grid_h {
            grid[gy * grid_w + gx] = Some(idx);
        }
        points.push((cx, cy));
        active.push(idx);
    }

    // Bridsonのメインループ
    while !active.is_empty() {
        // ランダムなアクティブ点を選択
        let active_idx = (next_rand() * active.len() as f64) as usize;
        let active_idx = active_idx.min(active.len() - 1);
        let pt_idx = active[active_idx];
        let (bx, by) = all_pts[pt_idx];

        let mut found = false;
        for _ in 0..max_attempts {
            // min_dist ～ 2*min_dist の環状領域にランダム点を生成
            let angle = next_rand() as f32 * std::f32::consts::TAU;
            let r = min_dist + next_rand() as f32 * min_dist;
            let nx = bx + angle.cos() * r;
            let ny = by + angle.sin() * r;

            if nx < 0.0 || nx >= w as f32 || ny < 0.0 || ny >= h as f32 {
                continue;
            }

            let ix = nx as usize;
            let iy = ny as usize;
            if ix >= w || iy >= h || alpha[iy * w + ix] < threshold {
                continue;
            }

            if is_valid_point(nx, ny, &all_pts, &grid, grid_w, grid_h, cell_size, min_dist2) {
                let new_idx = all_pts.len();
                all_pts.push((nx, ny));
                let gx = (nx / cell_size) as usize;
                let gy = (ny / cell_size) as usize;
                if gx < grid_w && gy < grid_h {
                    grid[gy * grid_w + gx] = Some(new_idx);
                }
                points.push((nx, ny));
                active.push(new_idx);
                found = true;
                break;
            }
        }

        if !found {
            active.swap_remove(active_idx);
        }
    }

    // 座標を画像中心原点に変換
    points.iter().map(|&(x, y)| (x - hw, y - hh)).collect()
}

fn is_valid_point(
    x: f32,
    y: f32,
    all_pts: &[(f32, f32)],
    grid: &[Option<usize>],
    grid_w: usize,
    grid_h: usize,
    cell_size: f32,
    min_dist2: f32,
) -> bool {
    let gx = (x / cell_size) as i32;
    let gy = (y / cell_size) as i32;

    // 近傍セルをチェック（5×5）
    for dy in -2..=2 {
        for dx in -2..=2 {
            let nx = gx + dx;
            let ny = gy + dy;
            if nx >= 0 && nx < grid_w as i32 && ny >= 0 && ny < grid_h as i32 {
                if let Some(idx) = grid[(ny as usize) * grid_w + (nx as usize)] {
                    let (px, py) = all_pts[idx];
                    let ddx = x - px;
                    let ddy = y - py;
                    if ddx * ddx + ddy * ddy < min_dist2 {
                        return false;
                    }
                }
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poisson_basic() {
        let w = 100;
        let h = 100;
        let alpha = vec![255u8; w * h];
        let min_dist = 10.0;
        let points = poisson_disk_sample(w, h, min_dist, &alpha, 128, &[], 42);
        assert!(!points.is_empty(), "点が生成されるべき");

        // すべての点間距離がmin_dist以上であることを確認
        let min_dist2 = min_dist * min_dist * 0.95; // 少しマージン
        for i in 0..points.len() {
            for j in (i + 1)..points.len() {
                let dx = points[i].0 - points[j].0;
                let dy = points[i].1 - points[j].1;
                let d2 = dx * dx + dy * dy;
                assert!(
                    d2 >= min_dist2,
                    "点{:?}と{:?}が近すぎる: dist={}",
                    points[i],
                    points[j],
                    d2.sqrt()
                );
            }
        }
    }

    #[test]
    fn test_poisson_alpha_region() {
        // 左半分のみ不透明
        let w = 100;
        let h = 100;
        let mut alpha = vec![0u8; w * h];
        for y in 0..h {
            for x in 0..w / 2 {
                alpha[y * w + x] = 255;
            }
        }
        let points = poisson_disk_sample(w, h, 10.0, &alpha, 128, &[], 42);
        let hw = w as f32 * 0.5;
        for &(px, _) in &points {
            // 中心原点なので左半分は px < 0
            assert!(
                px < 1.0,
                "右半分(透明領域)に点がある: x={}",
                px
            );
        }
    }
}
