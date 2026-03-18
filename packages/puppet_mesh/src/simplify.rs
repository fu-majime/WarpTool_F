/// Douglas-Peucker 輪郭単純化アルゴリズム
///
/// 輪郭線の点数を減らしながら、形状の特徴を保持する。
/// ※安全マージン（epsilon分の事前膨張）と組み合わせることで、
/// 　不透明領域を一切損なわず、かつピクセル階段を大幅にジャンプする
/// 　クリーンなポリゴンを生成する。

/// Douglas-Peucker アルゴリズムで輪郭を単純化する。
pub fn simplify(points: &[(f32, f32)], epsilon: f32) -> Vec<(f32, f32)> {
    let n = points.len();
    if n <= 2 {
        return points.to_vec();
    }

    let mut keep = vec![false; n];
    keep[0] = true;
    keep[n - 1] = true;

    dp_recursive(points, 0, n - 1, epsilon, &mut keep);

    points
        .iter()
        .enumerate()
        .filter(|(i, _)| keep[*i])
        .map(|(_, p)| *p)
        .collect()
}

/// 閉じた輪郭（ループ）を単純化する。
pub fn simplify_loop(points: &[(f32, f32)], epsilon: f32) -> Vec<(f32, f32)> {
    let n = points.len();
    if n <= 3 {
        return points.to_vec();
    }

    // 最も離れた2点を見つけて分割点にする
    let (split_a, split_b) = find_farthest_pair(points);

    let mut keep = vec![false; n];
    keep[split_a] = true;
    keep[split_b] = true;

    // 2つの区間に分けてDPを実行
    if split_a < split_b {
        dp_recursive(points, split_a, split_b, epsilon, &mut keep);
        dp_recursive_wrap(points, split_b, split_a, epsilon, &mut keep);
    } else {
        dp_recursive(points, split_b, split_a, epsilon, &mut keep);
        dp_recursive_wrap(points, split_a, split_b, epsilon, &mut keep);
    }

    points
        .iter()
        .enumerate()
        .filter(|(i, _)| keep[*i])
        .map(|(_, p)| *p)
        .collect()
}

fn find_farthest_pair(points: &[(f32, f32)]) -> (usize, usize) {
    let n = points.len();
    let mut max_d2 = 0.0f32;
    let mut pair = (0, 1);
    for i in 0..n {
        for j in (i + 1)..n {
            let d2 = dist2(points[i], points[j]);
            if d2 > max_d2 {
                max_d2 = d2;
                pair = (i, j);
            }
        }
    }
    pair
}

fn dp_recursive(
    points: &[(f32, f32)],
    start: usize,
    end: usize,
    epsilon: f32,
    keep: &mut [bool],
) {
    if end <= start + 1 {
        return;
    }

    let (max_idx, max_dist) = farthest_from_line(points, start, end);
    if max_dist > epsilon {
        keep[max_idx] = true;
        dp_recursive(points, start, max_idx, epsilon, keep);
        dp_recursive(points, max_idx, end, epsilon, keep);
    }
}

fn dp_recursive_wrap(
    points: &[(f32, f32)],
    start: usize,
    end: usize,
    epsilon: f32,
    keep: &mut [bool],
) {
    let n = points.len();
    let count = if end >= start {
        end - start + 1
    } else {
        (n - start) + end + 1
    };

    if count <= 2 {
        return;
    }

    let line_start = points[start];
    let line_end = points[end];
    let mut max_dist = 0.0f32;
    let mut max_idx = start;

    let mut i = (start + 1) % n;
    while i != end {
        let d = point_to_line_dist(points[i], line_start, line_end);
        if d > max_dist {
            max_dist = d;
            max_idx = i;
        }
        i = (i + 1) % n;
    }

    if max_dist > epsilon {
        keep[max_idx] = true;
        dp_recursive_wrap(points, start, max_idx, epsilon, keep);
        dp_recursive_wrap(points, max_idx, end, epsilon, keep);
    }
}

fn farthest_from_line(points: &[(f32, f32)], start: usize, end: usize) -> (usize, f32) {
    let a = points[start];
    let b = points[end];
    let mut max_dist = 0.0f32;
    let mut max_idx = start;

    for i in (start + 1)..end {
        let d = point_to_line_dist(points[i], a, b);
        if d > max_dist {
            max_dist = d;
            max_idx = i;
        }
    }

    (max_idx, max_dist)
}

fn point_to_line_dist(p: (f32, f32), a: (f32, f32), b: (f32, f32)) -> f32 {
    let dx = b.0 - a.0;
    let dy = b.1 - a.1;
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-10 {
        return ((p.0 - a.0).powi(2) + (p.1 - a.1).powi(2)).sqrt();
    }
    let cross = (p.0 - a.0) * dy - (p.1 - a.1) * dx;
    cross.abs() / len2.sqrt()
}

fn dist2(a: (f32, f32), b: (f32, f32)) -> f32 {
    let dx = a.0 - b.0;
    let dy = a.1 - b.1;
    dx * dx + dy * dy
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simplify_line() {
        let points: Vec<(f32, f32)> = (0..10).map(|i| (i as f32, 0.0)).collect();
        let result = simplify(&points, 0.1);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_simplify_preserves_corners() {
        let points = vec![
            (0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0),
            (4.0, 0.0), (5.0, 0.0), (5.0, 1.0), (5.0, 2.0),
            (5.0, 3.0), (5.0, 4.0), (5.0, 5.0),
        ];
        let result = simplify(&points, 0.1);
        assert_eq!(result.len(), 3);
    }
}
