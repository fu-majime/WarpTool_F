/// 制約付き Delaunay 三角形分割
///
/// Bowyer-Watson アルゴリズムで Delaunay 三角形分割を行い、
/// その後制約エッジ（輪郭辺）を挿入する。
/// 最後にα領域外の三角形を除去する。

use std::collections::HashSet;

const EPSILON: f64 = 1e-10;

#[derive(Clone)]
struct Triangle {
    v: [usize; 3], // 頂点インデックス
    cc: (f64, f64), // 外接円中心
    cr2: f64,       // 外接円半径の2乗
}

/// 制約付き Delaunay 三角形分割を実行する。
///
/// - `points`: 全頂点（輪郭点 + 内部点）、画像中心原点
/// - `contour_edges`: 輪郭の制約エッジ（頂点インデックスのペア）
/// - `alpha`: α画像データ
/// - `w`, `h`: 画像サイズ
/// - `threshold`: αしきい値
///
/// 戻り値: 三角形のリスト（各要素は3つの頂点インデックス、0-based）
pub fn triangulate(
    points: &[(f32, f32)],
    contour_edges: &[(usize, usize)],
    alpha: &[u8],
    w: usize,
    h: usize,
    threshold: u8,
) -> Vec<[usize; 3]> {
    let n = points.len();
    if n < 3 {
        return Vec::new();
    }

    let hw = w as f64 * 0.5;
    let hh = h as f64 * 0.5;

    // f64に変換
    let pts: Vec<(f64, f64)> = points.iter().map(|&(x, y)| (x as f64, y as f64)).collect();

    // バウンディングボックス
    let mut min_x = pts[0].0;
    let mut max_x = pts[0].0;
    let mut min_y = pts[0].1;
    let mut max_y = pts[0].1;
    for &(x, y) in &pts {
        if x < min_x { min_x = x; }
        if x > max_x { max_x = x; }
        if y < min_y { min_y = y; }
        if y > max_y { max_y = y; }
    }

    let dx = max_x - min_x;
    let dy = max_y - min_y;
    let dmax = dx.max(dy);
    let mid_x = (min_x + max_x) / 2.0;
    let mid_y = (min_y + max_y) / 2.0;

    // スーパー三角形: 全点を含む十分大きな三角形
    let ss = dmax * 20.0;
    let super_a = n;
    let super_b = n + 1;
    let super_c = n + 2;

    let mut all_pts = pts.clone();
    all_pts.push((mid_x - ss, mid_y - ss));
    all_pts.push((mid_x + ss, mid_y - ss));
    all_pts.push((mid_x, mid_y + ss));

    let cc = circumcircle(&all_pts, super_a, super_b, super_c);
    let mut triangles = vec![Triangle {
        v: [super_a, super_b, super_c],
        cc: (cc.0, cc.1),
        cr2: cc.2,
    }];

    // Bowyer-Watson: 各点を挿入
    for i in 0..n {
        let px = all_pts[i].0;
        let py = all_pts[i].1;

        let mut bad_tris: Vec<usize> = Vec::new();
        for (ti, tri) in triangles.iter().enumerate() {
            let ddx = px - tri.cc.0;
            let ddy = py - tri.cc.1;
            if ddx * ddx + ddy * ddy <= tri.cr2 + EPSILON {
                bad_tris.push(ti);
            }
        }

        // 多角形の穴のエッジを収集
        let mut polygon: Vec<(usize, usize)> = Vec::new();
        for &bi in &bad_tris {
            let tri = &triangles[bi];
            for edge_idx in 0..3 {
                let e = (tri.v[edge_idx], tri.v[(edge_idx + 1) % 3]);
                // 他の「悪い」三角形と共有されていないエッジのみ
                let shared = bad_tris.iter().any(|&oi| {
                    oi != bi && triangle_has_edge(&triangles[oi], e.0, e.1)
                });
                if !shared {
                    polygon.push(e);
                }
            }
        }

        // 悪い三角形を除去（逆順）
        bad_tris.sort_unstable();
        for &bi in bad_tris.iter().rev() {
            triangles.swap_remove(bi);
        }

        // 新しい三角形を作成
        for &(ea, eb) in &polygon {
            let cc = circumcircle(&all_pts, ea, eb, i);
            triangles.push(Triangle {
                v: [ea, eb, i],
                cc: (cc.0, cc.1),
                cr2: cc.2,
            });
        }
    }

    // スーパー三角形の頂点を含む三角形を除去
    triangles.retain(|tri| {
        tri.v[0] < n && tri.v[1] < n && tri.v[2] < n
    });

    // 制約エッジの挿入
    if !contour_edges.is_empty() {
        let _constraint_set: HashSet<(usize, usize)> = contour_edges
            .iter()
            .flat_map(|&(a, b)| vec![(a, b), (b, a)])
            .collect();

        // 制約エッジが存在しない場合、交差する三角形を再構成
        for &(ca, cb) in contour_edges {
            if ca >= n || cb >= n {
                continue;
            }
            // エッジが既に存在するかチェック
            let exists = triangles.iter().any(|tri| {
                triangle_has_edge(tri, ca, cb)
            });

            if !exists {
                // エッジフリップで制約エッジを挿入する
                insert_constraint_edge(&mut triangles, &all_pts, ca, cb);
            }
        }

        // α領域外の三角形を除去（重心チェック）
        triangles.retain(|tri| {
            let gx = (all_pts[tri.v[0]].0 + all_pts[tri.v[1]].0 + all_pts[tri.v[2]].0) / 3.0;
            let gy = (all_pts[tri.v[0]].1 + all_pts[tri.v[1]].1 + all_pts[tri.v[2]].1) / 3.0;
            let px = (gx + hw) as i32;
            let py = (gy + hh) as i32;
            if px >= 0 && px < w as i32 && py >= 0 && py < h as i32 {
                alpha[(py as usize) * w + (px as usize)] >= threshold
            } else {
                false
            }
        });
    } else {
        // 制約エッジなしの場合もα領域チェック
        triangles.retain(|tri| {
            let gx = (all_pts[tri.v[0]].0 + all_pts[tri.v[1]].0 + all_pts[tri.v[2]].0) / 3.0;
            let gy = (all_pts[tri.v[0]].1 + all_pts[tri.v[1]].1 + all_pts[tri.v[2]].1) / 3.0;
            let px = (gx + hw) as i32;
            let py = (gy + hh) as i32;
            if px >= 0 && px < w as i32 && py >= 0 && py < h as i32 {
                alpha[(py as usize) * w + (px as usize)] >= threshold
            } else {
                false
            }
        });
    }

    triangles.iter().map(|tri| tri.v).collect()
}

fn circumcircle(pts: &[(f64, f64)], a: usize, b: usize, c: usize) -> (f64, f64, f64) {
    let (ax, ay) = pts[a];
    let (bx, by) = pts[b];
    let (cx, cy) = pts[c];

    let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if d.abs() < EPSILON {
        let ux = (ax + bx + cx) / 3.0;
        let uy = (ay + by + cy) / 3.0;
        return (ux, uy, 1e18);
    }

    let a2 = ax * ax + ay * ay;
    let b2 = bx * bx + by * by;
    let c2 = cx * cx + cy * cy;

    let ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
    let uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;

    let r2 = (ax - ux) * (ax - ux) + (ay - uy) * (ay - uy);
    (ux, uy, r2)
}

fn triangle_has_edge(tri: &Triangle, a: usize, b: usize) -> bool {
    let v = &tri.v;
    for i in 0..3 {
        let va = v[i];
        let vb = v[(i + 1) % 3];
        if (va == a && vb == b) || (va == b && vb == a) {
            return true;
        }
    }
    false
}

/// 制約エッジを挿入する。
/// 交差する三角形のエッジをフリップして制約エッジを確立する。
fn insert_constraint_edge(
    triangles: &mut Vec<Triangle>,
    pts: &[(f64, f64)],
    ca: usize,
    cb: usize,
) {
    // 簡易実装: 制約エッジと交差するDelaunayエッジを見つけてフリップ
    // 複雑なケースでは失敗する可能性があるが、通常の輪郭には十分

    let max_iterations = 100;
    for _ in 0..max_iterations {
        // 現在の三角形群で制約エッジが存在するかチェック
        let exists = triangles.iter().any(|tri| triangle_has_edge(tri, ca, cb));
        if exists {
            return;
        }

        // 制約エッジと交差するDelaunayエッジを探す
        let mut flip_found = false;
        for i in 0..triangles.len() {
            for j in 0..3 {
                let ea = triangles[i].v[j];
                let eb = triangles[i].v[(j + 1) % 3];

                // 制約エッジの端点とDelaunayエッジの端点が同じ場合はスキップ
                if ea == ca || ea == cb || eb == ca || eb == cb {
                    continue;
                }

                if segments_intersect(pts[ca], pts[cb], pts[ea], pts[eb]) {
                    // 隣接三角形を探す
                    let mut adj_idx = None;
                    let mut adj_opp = 0;
                    for k in 0..triangles.len() {
                        if k == i {
                            continue;
                        }
                        if triangle_has_edge(&triangles[k], ea, eb) {
                            adj_idx = Some(k);
                            // 対向頂点を見つける
                            for m in 0..3 {
                                if triangles[k].v[m] != ea && triangles[k].v[m] != eb {
                                    adj_opp = triangles[k].v[m];
                                    break;
                                }
                            }
                            break;
                        }
                    }

                    if let Some(ai) = adj_idx {
                        // 対向頂点の取得
                        let opp_i = {
                            let mut o = 0;
                            for m in 0..3 {
                                if triangles[i].v[m] != ea && triangles[i].v[m] != eb {
                                    o = triangles[i].v[m];
                                    break;
                                }
                            }
                            o
                        };

                        // エッジフリップ: (ea, eb) → (opp_i, adj_opp)
                        let cc1 = circumcircle(pts, opp_i, adj_opp, ea);
                        let cc2 = circumcircle(pts, opp_i, adj_opp, eb);

                        // 大きいインデックスを先に除去
                        let (first_rm, second_rm) = if i > ai { (i, ai) } else { (ai, i) };
                        triangles.swap_remove(first_rm);
                        if second_rm < triangles.len() {
                            triangles.swap_remove(second_rm);
                        }

                        triangles.push(Triangle {
                            v: [opp_i, adj_opp, ea],
                            cc: (cc1.0, cc1.1),
                            cr2: cc1.2,
                        });
                        triangles.push(Triangle {
                            v: [opp_i, adj_opp, eb],
                            cc: (cc2.0, cc2.1),
                            cr2: cc2.2,
                        });

                        flip_found = true;
                        break;
                    }
                }
            }
            if flip_found {
                break;
            }
        }

        if !flip_found {
            break;
        }
    }
}

/// 2つの線分が交差するかを判定する（端点共有は除く）
fn segments_intersect(
    a1: (f64, f64),
    a2: (f64, f64),
    b1: (f64, f64),
    b2: (f64, f64),
) -> bool {
    let d1 = cross2d(b1, b2, a1);
    let d2 = cross2d(b1, b2, a2);
    let d3 = cross2d(a1, a2, b1);
    let d4 = cross2d(a1, a2, b2);

    if ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
        && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
    {
        return true;
    }

    false
}

fn cross2d(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_triangulate_square() {
        // 正方形の4頂点
        let points = vec![
            (-5.0f32, -5.0),
            (5.0, -5.0),
            (5.0, 5.0),
            (-5.0, 5.0),
        ];
        let w = 20;
        let h = 20;
        let alpha = vec![255u8; w * h];
        let tris = triangulate(&points, &[], &alpha, w, h, 1);
        assert_eq!(tris.len(), 2, "4点から2三角形が生成されるべき");
    }

    #[test]
    fn test_triangulate_with_interior() {
        // 5点: 正方形+中心
        let points = vec![
            (-10.0f32, -10.0),
            (10.0, -10.0),
            (10.0, 10.0),
            (-10.0, 10.0),
            (0.0, 0.0), // 中心
        ];
        let w = 30;
        let h = 30;
        let alpha = vec![255u8; w * h];
        let tris = triangulate(&points, &[], &alpha, w, h, 1);
        assert!(tris.len() >= 4, "5点から4三角形以上が生成されるべき");
    }
}
