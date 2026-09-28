//! Insert material handles into the rest topology, including both faces of a
//! shared edge. Exterior handles attach to the closest surface point; they
//! never add a long, transparent triangle extending out of the silhouette.
use aviutl2::anyhow::{self, Context as _};
type Point = (f64, f64);
fn cross(a: Point, b: Point, c: Point) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}
fn distance2(a: Point, b: Point) -> f64 {
    (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)
}

pub fn prepare_with_spacing(
    vertices: &[f64],
    indices: &[i32],
    pins: &[f64],
    radius: f64,
    spacing: f64,
) -> anyhow::Result<(Vec<f64>, Vec<i32>)> {
    anyhow::ensure!(
        spacing.is_finite() && spacing >= 0.0,
        "invalid refinement spacing"
    );
    type Cached = (Vec<f64>, Vec<i32>, Vec<f64>, f64, Vec<f64>, Vec<i32>, f64);
    thread_local! { static CACHE:std::cell::RefCell<std::collections::VecDeque<Cached>>=const {std::cell::RefCell::new(std::collections::VecDeque::new())}; }
    if let Some(result) = CACHE.with(|cache| {
        let mut cache = cache.borrow_mut();
        let i = cache.iter().position(|e| {
            e.0 == vertices && e.1 == indices && e.2 == pins && e.3 == radius && e.6 == spacing
        })?;
        let e = cache.remove(i)?;
        let result = (e.4.clone(), e.5.clone());
        cache.push_front(e);
        Some(result)
    }) {
        return Ok(result);
    }
    let result = prepare_uncached(vertices, indices, pins, radius, spacing)?;
    CACHE.with(|cache| {
        let mut cache = cache.borrow_mut();
        cache.push_front((
            vertices.to_vec(),
            indices.to_vec(),
            pins.to_vec(),
            radius,
            result.0.clone(),
            result.1.clone(),
            spacing,
        ));
        cache.truncate(4);
    });
    Ok(result)
}

fn prepare_uncached(
    vertices: &[f64],
    indices: &[i32],
    pins: &[f64],
    radius: f64,
    spacing: f64,
) -> anyhow::Result<(Vec<f64>, Vec<i32>)> {
    anyhow::ensure!(
        vertices.len().is_multiple_of(2)
            && pins.len().is_multiple_of(2)
            && indices.len().is_multiple_of(3),
        "invalid embedding arrays"
    );
    anyhow::ensure!(
        vertices.iter().chain(pins).all(|v| v.is_finite()) && radius.is_finite() && radius > 0.0,
        "invalid embedding coordinates"
    );
    let mut points = vertices
        .chunks_exact(2)
        .map(|v| (v[0], v[1]))
        .collect::<Vec<_>>();
    let mut triangles = indices
        .chunks_exact(3)
        .map(|t| {
            let a = usize::try_from(t[0]).context("negative index")?;
            let b = usize::try_from(t[1]).context("negative index")?;
            let c = usize::try_from(t[2]).context("negative index")?;
            anyhow::ensure!(
                [a, b, c].iter().all(|&i| i < points.len()),
                "invalid embedding index"
            );
            Ok([a, b, c])
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    if triangles.is_empty() {
        return Ok((vertices.to_vec(), indices.to_vec()));
    }
    let mut centers = Vec::new();
    let mut edge_counts = std::collections::BTreeMap::new();
    for t in &triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *edge_counts.entry((a.min(b), a.max(b))).or_insert(0_usize) += 1;
        }
    }
    let boundary_segments = edge_counts
        .into_iter()
        .filter_map(|((a, b), n)| (n == 1).then_some((points[a], points[b])))
        .collect::<Vec<_>>();
    for pin in pins.chunks_exact(2) {
        let p = (pin[0], pin[1]);
        let p = if contains(&points, &triangles, p) {
            p
        } else {
            project(&points, &triangles, p)
        };
        insert(&mut points, &mut triangles, p);
        let thickness = boundary_segments
            .iter()
            .map(|&(a, b)| {
                let (x, y) = (b.0 - a.0, b.1 - a.1);
                let t = (((p.0 - a.0) * x + (p.1 - a.1) * y) / (x * x + y * y).max(1e-20))
                    .clamp(0.0, 1.0);
                (p.0 - a.0 - t * x).hypot(p.1 - a.1 - t * y)
            })
            .fold(f64::INFINITY, f64::min);
        centers.push((p, thickness.max(radius).min(radius * 8.0) * 0.5));
    }
    // Local refinement provides a finite neighbourhood for position and detail
    // handles without increasing density across the whole object.
    for &(p, radius) in &centers {
        for k in 0..4 {
            let angle = k as f64 * std::f64::consts::TAU / 4.0;
            let q = (p.0 + radius * angle.cos(), p.1 + radius * angle.sin());
            if contains(&points, &triangles, q)
                && boundary_segments.iter().all(|&(a, b)| {
                    let length2 = distance2(a, b);
                    let t = (((q.0 - a.0) * (b.0 - a.0) + (q.1 - a.1) * (b.1 - a.1))
                        / length2.max(1e-20))
                    .clamp(0., 1.);
                    distance2(q, (a.0 + t * (b.0 - a.0), a.1 + t * (b.1 - a.1)))
                        > (radius * 0.2).max(spacing * 0.1).powi(2)
                })
                && points
                    .iter()
                    .all(|&other| distance2(q, other) > (radius * 0.35).max(spacing * 0.15).powi(2))
            {
                insert(&mut points, &mut triangles, q);
            }
        }
    }
    // Restore Delaunay interior edges after insertion instead of keeping thin
    // triangle fans. Only the material boundary is a triangulation constraint.
    let mut counts = std::collections::BTreeMap::<(usize, usize), usize>::new();
    for t in &triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *counts.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let mut boundary = counts
        .into_iter()
        .filter_map(|(edge, n)| (n == 1).then_some(edge))
        .collect::<Vec<_>>();
    // A mandatory pin may land almost on an existing interior sample. Keep
    // the exact pin and the constrained outline, then retriangulate without
    // the redundant sample rather than making an almost-zero edge.
    let outline = boundary
        .iter()
        .flat_map(|&(a, b)| [a, b])
        .collect::<std::collections::BTreeSet<_>>();
    let clearance2 = (spacing * 0.15).max(radius * 0.2).powi(2);
    let mut remap = vec![usize::MAX; points.len()];
    let mut retained = Vec::new();
    for (i, &p) in points.iter().enumerate() {
        let is_pin = centers
            .iter()
            .any(|&(center, _)| distance2(p, center) < 1e-16);
        let redundant = !outline.contains(&i)
            && !is_pin
            && centers
                .iter()
                .any(|&(center, _)| distance2(p, center) < clearance2);
        if !redundant {
            remap[i] = retained.len();
            retained.push(p);
        }
    }
    for (a, b) in &mut boundary {
        *a = remap[*a];
        *b = remap[*b];
    }
    points = retained;
    let triangles = crate::delaunay::triangulate_f64(&points, &boundary)?;
    Ok((
        points.into_iter().flat_map(|p| [p.0, p.1]).collect(),
        triangles.into_iter().flatten().map(|i| i as i32).collect(),
    ))
}
fn bary(points: &[Point], t: [usize; 3], p: Point) -> Option<[f64; 3]> {
    let [a, b, c] = t.map(|i| points[i]);
    let area = cross(a, b, c);
    if area.abs() < 1.0e-14 {
        return None;
    }
    let b = [
        cross(b, c, p) / area,
        cross(c, a, p) / area,
        cross(a, b, p) / area,
    ];
    b.iter().all(|&v| v >= -1.0e-10).then_some(b)
}
fn contains(points: &[Point], triangles: &[[usize; 3]], p: Point) -> bool {
    triangles.iter().any(|&t| bary(points, t, p).is_some())
}
fn project(points: &[Point], triangles: &[[usize; 3]], p: Point) -> Point {
    let mut best = points[triangles[0][0]];
    let mut d = distance2(p, best);
    for t in triangles {
        for k in 0..3 {
            let a = points[t[k]];
            let b = points[t[(k + 1) % 3]];
            let length = distance2(a, b);
            if length == 0.0 {
                continue;
            }
            let s =
                (((p.0 - a.0) * (b.0 - a.0) + (p.1 - a.1) * (b.1 - a.1)) / length).clamp(0.0, 1.0);
            let q = (a.0 + s * (b.0 - a.0), a.1 + s * (b.1 - a.1));
            let candidate = distance2(p, q);
            if candidate < d {
                d = candidate;
                best = q;
            }
        }
    }
    best
}
fn insert(points: &mut Vec<Point>, triangles: &mut Vec<[usize; 3]>, p: Point) {
    if points.iter().any(|&q| distance2(p, q) < 1.0e-16) {
        return;
    }
    let Some((index, b)) = triangles
        .iter()
        .enumerate()
        .find_map(|(i, &t)| bary(points, t, p).map(|b| (i, b)))
    else {
        return;
    };
    let t = triangles[index];
    let n = points.len();
    points.push(p);
    if let Some(k) = b.iter().position(|v| v.abs() < 1.0e-10) {
        let (a, b) = (t[(k + 1) % 3], t[(k + 2) % 3]);
        let mut updated = Vec::with_capacity(triangles.len() + 2);
        for &t in triangles.iter() {
            if t.contains(&a) && t.contains(&b) {
                for j in 0..3 {
                    if (t[j] == a && t[(j + 1) % 3] == b) || (t[j] == b && t[(j + 1) % 3] == a) {
                        updated.push([t[j], n, t[(j + 2) % 3]]);
                        updated.push([n, t[(j + 1) % 3], t[(j + 2) % 3]]);
                        break;
                    }
                }
            } else {
                updated.push(t);
            }
        }
        *triangles = updated;
    } else {
        triangles[index] = [t[0], t[1], n];
        triangles.push([t[1], t[2], n]);
        triangles.push([t[2], t[0], n]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_pin_replaces_a_nearly_coincident_interior_sample() {
        let (v, t) = prepare_with_spacing(
            &[0., 0., 100., 0., 100., 100., 0., 100., 50.0001, 50.],
            &[0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4],
            &[50., 50.],
            8.,
            30.,
        )
        .unwrap();
        assert!(v.chunks_exact(2).any(|p| p == [50., 50.]));
        assert!(!v.chunks_exact(2).any(|p| p == [50.0001, 50.]));
        let area: f64 = t
            .chunks_exact(3)
            .map(|t| {
                let p = |i: i32| (v[i as usize * 2], v[i as usize * 2 + 1]);
                cross(p(t[0]), p(t[1]), p(t[2])) * 0.5
            })
            .sum();
        assert!((area - 10000.).abs() < 1e-8);
    }
    #[test]
    fn optional_handle_samples_keep_clear_of_long_boundary_edges() {
        let (v, t) = prepare_with_spacing(
            &[0., 0., 100., 0., 100., 100., 0., 100.],
            &[0, 1, 2, 0, 2, 3],
            &[50., 5.],
            8.,
            30.,
        )
        .unwrap();
        assert!(v.chunks_exact(2).any(|p| p == [50., 5.]));
        for p in v.chunks_exact(2).skip(4) {
            assert!(
                p[0].min(p[1]).min(100. - p[0]).min(100. - p[1]) >= 3.,
                "optional sample creates a boundary sliver: {p:?}"
            );
        }
        let area: f64 = t
            .chunks_exact(3)
            .map(|t| {
                let p = |i: i32| (v[i as usize * 2], v[i as usize * 2 + 1]);
                cross(p(t[0]), p(t[1]), p(t[2])) * 0.5
            })
            .sum();
        assert!((area - 10000.).abs() < 1e-8);
    }
    #[test]
    fn shared_edge_is_split_without_cracks_or_area_loss() {
        let (v, t) = prepare(
            &[0., 0., 10., 0., 10., 10., 0., 10.],
            &[0, 1, 2, 0, 2, 3],
            &[5., 5.],
            2.,
        )
        .unwrap();
        let area = t
            .chunks_exact(3)
            .map(|t| {
                let p = |i: i32| (v[i as usize * 2], v[i as usize * 2 + 1]);
                let area = cross(p(t[0]), p(t[1]), p(t[2])) * 0.5;
                assert!(area > 1.0e-10);
                area
            })
            .sum::<f64>();
        assert!((area - 100.0).abs() < 1.0e-8);
        assert!(v.chunks_exact(2).any(|p| p == [5., 5.]));
    }
}
