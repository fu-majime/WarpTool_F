/// 高速縁取り膨張（Euclidean Distance Transformベース）
///
/// α二値画像に対して、外側方向にのみ正確に膨張させる。
/// Meijster/Roerdink/Hesselink の O(w×h) EDT アルゴリズムを使用。

/// 二値マップを指定半径だけ外側に膨張する。
///
/// - `binary`: 入力二値マップ (true = 不透明), サイズ w*h, **in-place で書き換えられる**
/// - `w`, `h`: 画像サイズ
/// - `radius`: 膨張半径（ピクセル）
///
/// O(w×h) の 2パス Euclidean Distance Transform によって
/// 各透明ピクセルから最近接不透明ピクセルへの2乗距離を計算し、
/// radius² 以下の透明ピクセルを不透明に変換する。
/// 元から不透明だったピクセルは一切変更しない（内側に侵食しない）。
pub fn dilate_edt(binary: &mut [bool], w: usize, h: usize, radius: f32) {
    if radius <= 0.0 || w == 0 || h == 0 {
        return;
    }

    let r2 = (radius * radius) as u64;
    let inf: u32 = (w + h) as u32; // 十分大きな値（1D距離として）

    // パス1: 各列について、最近接の不透明ピクセルまでの距離(の2乗)を1D計算
    // g[y * w + x] = 列x方向の最近接不透明ピクセルまでの距離
    let mut g = vec![0u32; w * h];

    // 横方向パス: 各行について
    for y in 0..h {
        // 左→右: 最近接不透明ピクセルまでの距離
        let mut dt = inf;
        for x in 0..w {
            if binary[y * w + x] {
                dt = 0;
            } else {
                dt = dt.saturating_add(1);
            }
            g[y * w + x] = dt;
        }
        // 右→左
        dt = inf;
        for x in (0..w).rev() {
            if binary[y * w + x] {
                dt = 0;
            } else {
                dt = dt.saturating_add(1);
            }
            g[y * w + x] = g[y * w + x].min(dt);
        }
    }

    // パス2: 各行について、g値を使って2D Euclidean距離の2乗を計算
    let mut s = vec![0i32; h]; // 包絡線の頂点位置
    let mut t = vec![0i32; h]; // 包絡線の列位置
    let mut result = vec![false; w * h];

    // 元の不透明ピクセルをコピー
    result.copy_from_slice(binary);

    for x in 0..w {
        // Meijster の列方向パス
        let mut q: i32 = 0; // 包絡線のスタック頂点数 - 1

        s[0] = 0;
        t[0] = 0;

        // f(x, i) = g[i*w+x]² + (y-i)²  の下側包絡線を構築
        let f = |i: i32, y: i32| -> u64 {
            let gi = g[i as usize * w + x] as u64;
            let dy = (y - i).unsigned_abs() as u64;
            gi * gi + dy * dy
        };

        // Sep(i, u): f(x,i) と f(x,u) の交点
        let sep = |i: i32, u: i32| -> i32 {
            let gi2 = (g[i as usize * w + x] as i64) * (g[i as usize * w + x] as i64);
            let gu2 = (g[u as usize * w + x] as i64) * (g[u as usize * w + x] as i64);
            let num = gu2 - gi2 + (u as i64) * (u as i64) - (i as i64) * (i as i64);
            let den = 2 * (u as i64 - i as i64);
            if den == 0 {
                return i32::MAX;
            }
            // 切り上げ除算
            (num / den) as i32
        };

        for u in 1..h as i32 {
            while q >= 0 && f(s[q as usize], t[q as usize]) >= f(u, t[q as usize]) {
                q -= 1;
            }
            if q < 0 {
                q = 0;
                s[0] = u;
            } else {
                let w_val = sep(s[q as usize], u) + 1;
                if w_val < h as i32 {
                    q += 1;
                    s[q as usize] = u;
                    t[q as usize] = w_val;
                }
            }
        }

        // 包絡線をスキャンして距離を判定
        for y in (0..h as i32).rev() {
            let dist2 = f(s[q as usize], y);
            if dist2 <= r2 {
                result[y as usize * w + x] = true;
            }
            if y == t[q as usize] && q > 0 {
                q -= 1;
            }
        }
    }

    binary.copy_from_slice(&result);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_dilation() {
        let mut binary = vec![false; 10 * 10];
        binary[55] = true; // (5,5)
        let original = binary.clone();
        dilate_edt(&mut binary, 10, 10, 0.0);
        assert_eq!(binary, original, "半径0では変化なし");
    }

    #[test]
    fn test_single_pixel_dilation() {
        let w = 20;
        let h = 20;
        let mut binary = vec![false; w * h];
        binary[10 * w + 10] = true; // 中央 (10,10)
        dilate_edt(&mut binary, w, h, 3.0);

        // (10,10) から半径3以内のピクセルが true になっているか
        for y in 0..h {
            for x in 0..w {
                let dx = x as f32 - 10.0;
                let dy = y as f32 - 10.0;
                let dist = (dx * dx + dy * dy).sqrt();
                if dist <= 3.0 {
                    assert!(binary[y * w + x], "({},{}) は膨張範囲内だが false", x, y);
                }
            }
        }
    }

    #[test]
    fn test_preserves_opaque() {
        let w = 10;
        let h = 10;
        let mut binary = vec![true; w * h];
        dilate_edt(&mut binary, w, h, 5.0);
        assert!(binary.iter().all(|&b| b), "全面不透明は変化なし");
    }
}
