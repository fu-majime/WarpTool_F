/// Marching Squares 輪郭抽出
///
/// α画像からしきい値に基づいて等値線(輪郭)を抽出する。
/// 結果は連結された輪郭ループのリストとして返される。

/// α値の2Dグリッド上でMarching Squaresを実行し、輪郭ループを抽出する。
///
/// - `alpha`: 各ピクセルのα値 (0-255), サイズ w*h
/// - `w`, `h`: 画像サイズ
/// - `threshold`: α境界しきい値 (0-255)
/// - `dilate`: 輪郭を外側に拡張するピクセル数
///
/// 戻り値: 輪郭ループのリスト。各ループは (x, y) 座標の列（画像中心原点）。
pub fn extract_contours(
    alpha: &[u8],
    w: usize,
    h: usize,
    threshold: u8,
    dilate: f32,
) -> Vec<Vec<(f32, f32)>> {
    if w < 2 || h < 2 {
        return Vec::new();
    }

    let hw = w as f32 * 0.5;
    let hh = h as f32 * 0.5;

    // 二値化マップ: true = 不透明領域
    let mut binary = vec![false; w * h];
    for i in 0..(w * h) {
        binary[i] = alpha[i] >= threshold;
    }

    // 膨張処理: 不透明領域を dilate ピクセル分拡張して
    // 輪郭単純化後にエッジが欠けるのを防ぐ
    // 高速化: 不透明ピクセルからのスタンプ方式 O(opaque_count × r²)
    if dilate > 0.0 {
        let radius = dilate.ceil() as i32;
        let r2 = dilate * dilate;
        let mut expanded = binary.clone();

        // 不透明ピクセルの近傍をマーク（逆方向スタンプ）
        for y in 0..h as i32 {
            for x in 0..w as i32 {
                if !binary[(y as usize) * w + (x as usize)] {
                    continue; // 透明ピクセルはスキップ
                }
                // この不透明ピクセルの周囲をマーク
                let y_start = (y - radius).max(0);
                let y_end = (y + radius).min(h as i32 - 1);
                let x_start = (x - radius).max(0);
                let x_end = (x + radius).min(w as i32 - 1);
                for ny in y_start..=y_end {
                    for nx in x_start..=x_end {
                        let dx = nx - x;
                        let dy = ny - y;
                        if (dx * dx + dy * dy) as f32 <= r2 {
                            expanded[(ny as usize) * w + (nx as usize)] = true;
                        }
                    }
                }
            }
        }
        binary = expanded;
    }
    marching_squares(&binary, w, h)
}

/// 事前に二値化済みのマップから Marching Squares で輪郭を抽出する。
/// EDT膨張との組み合わせで使用。
pub fn extract_contours_from_binary(
    binary: &[bool],
    w: usize,
    h: usize,
) -> Vec<Vec<(f32, f32)>> {
    if w < 2 || h < 2 {
        return Vec::new();
    }
    marching_squares(binary, w, h)
}

/// Marching Squares 実行（内部関数）
fn marching_squares(binary: &[bool], w: usize, h: usize) -> Vec<Vec<(f32, f32)>> {
    let hw = w as f32 * 0.5;
    let hh = h as f32 * 0.5;

    // Marching Squares: セル(x, y)の4隅の状態から等値線セグメントを抽出
    // セルは (w-1)×(h-1) 個
    // 各セグメントは (x0,y0) → (x1,y1) のエッジ
    let mut segments: Vec<((f32, f32), (f32, f32))> = Vec::new();

    for cy in 0..(h - 1) {
        for cx in 0..(w - 1) {
            // 4隅の値 (左上, 右上, 右下, 左下)
            let tl = binary[cy * w + cx];
            let tr = binary[cy * w + cx + 1];
            let br = binary[(cy + 1) * w + cx + 1];
            let bl = binary[(cy + 1) * w + cx];

            let case = (tl as u8) << 3 | (tr as u8) << 2 | (br as u8) << 1 | (bl as u8);

            if case == 0 || case == 15 {
                continue; // 全部同じ → エッジなし
            }

            // セル座標系: (cx, cy) が左上, (cx+1, cy+1) が右下
            // エッジ中点の座標 (線形補間は不要: 二値化なので中点でOK)
            let top = (cx as f32 + 0.5, cy as f32);
            let right = (cx as f32 + 1.0, cy as f32 + 0.5);
            let bottom = (cx as f32 + 0.5, cy as f32 + 1.0);
            let left = (cx as f32, cy as f32 + 0.5);

            // 16パターンのルックアップ
            match case {
                1 => segments.push((bottom, left)),
                2 => segments.push((right, bottom)),
                3 => segments.push((right, left)),
                4 => segments.push((top, right)),
                5 => {
                    // サドルポイント: 交差
                    segments.push((top, left));
                    segments.push((right, bottom));
                }
                6 => segments.push((top, bottom)),
                7 => segments.push((top, left)),
                8 => segments.push((left, top)),
                9 => segments.push((bottom, top)),
                10 => {
                    // サドルポイント: 交差
                    segments.push((left, bottom));
                    segments.push((top, right));
                }
                11 => segments.push((right, top)),
                12 => segments.push((left, right)),
                13 => segments.push((bottom, right)),
                14 => segments.push((left, bottom)),
                _ => {}
            }
        }
    }

    // セグメントを連結して輪郭ループにする
    let contours = chain_segments(&segments, hw, hh);
    contours
}

/// セグメントリストを連結してループ（または開いた輪郭線）にする。
/// 座標は画像中心原点に変換される。
fn chain_segments(
    segments: &[((f32, f32), (f32, f32))],
    hw: f32,
    hh: f32,
) -> Vec<Vec<(f32, f32)>> {
    let n = segments.len();
    if n == 0 {
        return Vec::new();
    }

    let mut used = vec![false; n];
    let mut contours = Vec::new();

    // 終点→セグメントインデックスのマルチマップ
    // 近接点をスナップするために量子化キーを使う
    let quantize = |x: f32, y: f32| -> (i32, i32) {
        ((x * 4.0).round() as i32, (y * 4.0).round() as i32)
    };

    // 始点から探索可能なセグメントのインデックスマップ
    let mut start_map: std::collections::HashMap<(i32, i32), Vec<usize>> =
        std::collections::HashMap::new();
    for (i, seg) in segments.iter().enumerate() {
        let key = quantize(seg.0 .0, seg.0 .1);
        start_map.entry(key).or_default().push(i);
    }

    for seed in 0..n {
        if used[seed] {
            continue;
        }
        used[seed] = true;

        let mut chain: Vec<(f32, f32)> = Vec::new();
        chain.push((segments[seed].0 .0 - hw, segments[seed].0 .1 - hh));
        chain.push((segments[seed].1 .0 - hw, segments[seed].1 .1 - hh));

        // 末尾から連結
        loop {
            let last = chain.last().unwrap();
            let key = quantize(last.0 + hw, last.1 + hh);
            let mut found = false;
            if let Some(indices) = start_map.get(&key) {
                for &idx in indices {
                    if !used[idx] {
                        used[idx] = true;
                        chain.push((segments[idx].1 .0 - hw, segments[idx].1 .1 - hh));
                        found = true;
                        break;
                    }
                }
            }
            if !found {
                break;
            }
        }

        if chain.len() >= 3 {
            contours.push(chain);
        }
    }

    contours
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_square() {
        // 10x10 の画像、中央4x4が不透明
        let w = 10;
        let h = 10;
        let mut alpha = vec![0u8; w * h];
        for y in 3..7 {
            for x in 3..7 {
                alpha[y * w + x] = 255;
            }
        }
        let contours = extract_contours(&alpha, w, h, 128, 0.0);
        assert!(!contours.is_empty(), "正方形から輪郭が抽出されるべき");
        // 輪郭点数は最低4点以上
        let total_pts: usize = contours.iter().map(|c| c.len()).sum();
        assert!(total_pts >= 4, "輪郭点が少なすぎる: {}", total_pts);
    }

    #[test]
    fn test_empty_image() {
        let w = 10;
        let h = 10;
        let alpha = vec![0u8; w * h];
        let contours = extract_contours(&alpha, w, h, 128, 0.0);
        assert!(contours.is_empty(), "空画像では輪郭なし");
    }

    #[test]
    fn test_full_image() {
        let w = 10;
        let h = 10;
        let alpha = vec![255u8; w * h];
        let contours = extract_contours(&alpha, w, h, 128, 0.0);
        assert!(contours.is_empty(), "全面不透明では内部輪郭なし");
    }

    #[test]
    fn test_dilation() {
        let w = 10;
        let h = 10;
        let mut alpha = vec![0u8; w * h];
        // 単一ピクセル
        alpha[5 * w + 5] = 255;
        let no_dilate = extract_contours(&alpha, w, h, 128, 0.0);
        let with_dilate = extract_contours(&alpha, w, h, 128, 2.0);
        // 膨張あり版のほうが輪郭点が多い（広い領域をカバー）
        let pts_no: usize = no_dilate.iter().map(|c| c.len()).sum();
        let pts_yes: usize = with_dilate.iter().map(|c| c.len()).sum();
        assert!(
            pts_yes >= pts_no,
            "膨張版の輪郭点数({})が非膨張版({})より少ない",
            pts_yes,
            pts_no
        );
    }
}
