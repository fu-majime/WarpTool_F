//! Density-scaled contour approximation. RDP chooses spans; each span is moved
//! outward only as far as its own omitted arc requires. Unlike one-sided vertex
//! deletion, diagonal pixel staircases do not survive as hundreds of handles.
type P = (f64, f64);

#[cfg(test)]
#[test]
fn only_subresolution_open_notches_are_bridged() {
    let outline = |depth| {
        vec![
            (0., 0.),
            (40., 0.),
            (40., 40.),
            (24., 40.),
            (20., 40. - depth),
            (16., 40.),
            (0., 40.),
        ]
    };
    let shallow = outline(6.);
    let cleaned = clean_short_steps(&shallow, 1., 4., 100., &[]);
    assert!(cleaned.len() < shallow.len());
    assert!(valid(&shallow, &cleaned));
    let deep = outline(20.);
    assert_eq!(clean_short_steps(&deep, 1., 4., 100., &[]), deep);
    assert_eq!(
        clean_short_steps(&shallow, 1., 4., 100., &[(20., 38.)]),
        shallow
    );
    let hole = shallow.iter().copied().rev().collect::<Vec<_>>();
    assert_eq!(clean_short_steps(&hole, 1., 4., 100., &[]), hole);
}

#[cfg(test)]
#[test]
fn short_mixed_turn_step_is_removed_without_cutting_material() {
    let raw = [
        (0., 0.),
        (100., 0.),
        (100., 50.),
        (99., 50.),
        (99., 100.),
        (0., 100.),
    ];
    let cleaned = clean_short_steps(&raw, 4., 4., 60., &[]);
    assert!(cleaned.len() < raw.len());
    assert!(valid(&raw, &cleaned));
    for i in 0..cleaned.len() {
        assert!(dist2(cleaned[i], cleaned[(i + 1) % cleaned.len()]) >= 16.);
    }
}

// Mixed convex/concave raster steps cannot be contracted by the convex-only
// supporting-line rule. Try bounded local replacements and explicitly verify
// containment and topology before accepting one.
fn clean_short_steps(
    raw: &[P],
    minimum: f64,
    tolerance: f64,
    budget: f64,
    blockers: &[P],
) -> Vec<P> {
    let mut polygon = raw.to_vec();
    let original_area = area(raw);
    loop {
        let n = polygon.len();
        if n <= 3 {
            break;
        }
        let mut best: Option<(f64, Vec<P>)> = None;
        for i in 0..n {
            let j = (i + 1) % n;
            let (a, b, c, d) = (
                polygon[(i + n - 1) % n],
                polygon[i],
                polygon[j],
                polygon[(j + 1) % n],
            );
            let mut candidates = Vec::new();
            let mouth2 = dist2(a, c);
            let turn = cross(a, b, c);
            // A narrow shallow indentation below the contour resolution does
            // not need two facing boundary fans. Bridge it in the calculation
            // mesh only; texture alpha still describes the visible gap.
            // Deep cuts and closed holes are deliberately excluded.
            if original_area > 0.
                && turn < 0.
                && mouth2 > 1e-12
                && mouth2 <= 16. * tolerance * tolerance
                && turn * turn / mouth2 <= 4. * tolerance * tolerance
            {
                candidates.push(
                    (0..n)
                        .filter(|&k| k != i)
                        .map(|k| polygon[k])
                        .collect::<Vec<_>>(),
                );
            }
            if dist2(b, c) < minimum * minimum {
                let mut replacements = vec![b, c];
                let u = (b.0 - a.0, b.1 - a.1);
                let v = (d.0 - c.0, d.1 - c.1);
                let det = u.0 * v.1 - u.1 * v.0;
                if det.abs() > 1e-12 {
                    let t = ((c.0 - a.0) * v.1 - (c.1 - a.1) * v.0) / det;
                    replacements.push((a.0 + t * u.0, a.1 + t * u.1));
                }
                for q in replacements {
                    let q = (q.0 as f32 as f64, q.1 as f32 as f64);
                    if dist2(q, b).max(dist2(q, c)) <= tolerance * tolerance {
                        candidates.push(
                            (0..n)
                                .filter(|&k| k != j)
                                .map(|k| if k == i { q } else { polygon[k] })
                                .collect(),
                        );
                    }
                }
            }
            for candidate in candidates {
                // Reject candidates outside the budget or worse than the best
                // before the quadratic duplicate/topology checks.
                let gain = (area(&candidate) - original_area).abs();
                if gain > budget || best.as_ref().is_some_and(|(cost, _)| gain >= *cost) {
                    continue;
                }
                if candidate
                    .iter()
                    .enumerate()
                    .any(|(k, &p)| candidate[k + 1..].iter().any(|&r| dist2(p, r) < 1e-10))
                {
                    continue;
                }
                if valid(&polygon, &candidate)
                    && !blockers
                        .iter()
                        .any(|&p| inside(p, &candidate) && !inside(p, &polygon))
                {
                    best = Some((gain, candidate));
                }
            }
        }
        match best {
            Some((_, next)) => polygon = next,
            None => break,
        }
    }
    polygon
}
// Monotone outward edge contraction. Each accepted operation adds material;
// no whole-contour offset and no inward chord across opaque pixels is needed.
fn contract(raw: &[P], tolerance: f64, budget: f64, blockers: &[P]) -> Vec<P> {
    contract_short_edges(raw, tolerance, budget, blockers, f64::INFINITY)
}
fn contract_short_edges(
    raw: &[P],
    tolerance: f64,
    budget: f64,
    blockers: &[P],
    short_edge: f64,
) -> Vec<P> {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    let n = raw.len();
    let mut p = raw.to_vec();
    let mut prev = (0..n).map(|i| (i + n - 1) % n).collect::<Vec<_>>();
    let mut next = (0..n).map(|i| (i + 1) % n).collect::<Vec<_>>();
    let mut alive = vec![true; n];
    let mut versions = vec![0usize; n];
    let mut sorted = raw.to_vec();
    sorted.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    let mut hull = Vec::new();
    for &q in &sorted {
        while hull.len() >= 2 && cross(hull[hull.len() - 2], hull[hull.len() - 1], q) <= 0. {
            hull.pop();
        }
        hull.push(q);
    }
    let lower = hull.len();
    for &q in sorted.iter().rev().skip(1) {
        while hull.len() > lower && cross(hull[hull.len() - 2], hull[hull.len() - 1], q) <= 0. {
            hull.pop();
        }
        hull.push(q);
    }
    hull.pop();
    let proposal = |i: usize, p: &[P], prev: &[usize], next: &[usize]| -> Option<(f64, P, bool)> {
        let a = p[prev[i]];
        let b = p[i];
        let c = p[next[i]];
        let turn = cross(a, b, c);
        if turn <= 1e-10 {
            if dist2(a, b).min(dist2(b, c)) > short_edge * short_edge {
                return None;
            }
            if short_edge.is_finite() {
                let depth = (0..hull.len())
                    .map(|j| {
                        let u = hull[j];
                        let v = hull[(j + 1) % hull.len()];
                        let t = (((b.0 - u.0) * (v.0 - u.0) + (b.1 - u.1) * (v.1 - u.1))
                            / dist2(u, v).max(1e-12))
                        .clamp(0., 1.);
                        dist2(b, (u.0 + t * (v.0 - u.0), u.1 + t * (v.1 - u.1)))
                    })
                    .fold(f64::INFINITY, f64::min);
                if depth > tolerance * tolerance {
                    return None;
                }
            }
            return Some((-turn, c, false));
        }
        if dist2(b, c) > short_edge * short_edge {
            return None;
        }
        let d = p[next[next[i]]];
        if cross(b, c, d) <= 0.0 {
            return None;
        }
        let u = (b.0 - a.0, b.1 - a.1);
        let v = (d.0 - c.0, d.1 - c.1);
        let det = u.0 * v.1 - u.1 * v.0;
        if det.abs() < 1e-10 {
            return None;
        }
        let t = ((c.0 - a.0) * v.1 - (c.1 - a.1) * v.0) / det;
        let q = (a.0 + t * u.0, a.1 + t * u.1);
        if t < 1.0 || dist2(q, b).max(dist2(q, c)) > tolerance * tolerance {
            return None;
        }
        let gain = cross(b, q, c);
        (gain >= -1e-8).then_some((gain.max(0.0), q, true))
    };
    let mut heap = BinaryHeap::new();
    for i in 0..n {
        if let Some((cost, _, _)) = proposal(i, &p, &prev, &next) {
            heap.push(Reverse((cost.to_bits(), i, 0usize)));
        }
    }
    let mut count = n;
    let mut spent = 0.0;
    while let Some(Reverse((_, i, version))) = heap.pop() {
        if !alive[i] || version != versions[i] || count <= 3 {
            continue;
        }
        let Some((cost, q, merge)) = proposal(i, &p, &prev, &next) else {
            continue;
        };
        if spent + cost > budget || cost > tolerance * tolerance {
            continue;
        }
        let a = prev[i];
        let b = next[i];
        let c = next[b];
        let added = if merge {
            [p[i], q, p[b]]
        } else {
            [p[a], p[b], p[i]]
        };
        if blockers.iter().any(|&point| inside(point, &added)) {
            continue;
        }
        let edges = if merge {
            vec![(p[a], q), (q, p[c])]
        } else {
            vec![(p[a], p[b])]
        };
        let intersects = (0..n)
            .filter(|&j| alive[j] && j != a && j != i && (!merge || j != b))
            .any(|j| {
                edges
                    .iter()
                    .any(|&(u, v)| proper_cross(u, v, p[j], p[next[j]]))
            });
        if intersects {
            continue;
        }
        if merge {
            p[i] = q;
            next[i] = c;
            prev[c] = i;
            alive[b] = false;
        } else {
            next[a] = b;
            prev[b] = a;
            alive[i] = false;
        }
        count -= 1;
        spent += cost;
        let center = if merge { i } else { a };
        for j in [prev[prev[center]], prev[center], center, next[center]] {
            versions[j] += 1;
            if let Some((cost, _, _)) = proposal(j, &p, &prev, &next) {
                heap.push(Reverse((cost.to_bits(), j, versions[j])));
            }
        }
    }
    let first = alive.iter().position(|&v| v).unwrap();
    let mut result = vec![p[first]];
    let mut i = next[first];
    while i != first {
        result.push(p[i]);
        i = next[i];
    }
    result
}
fn cross(a: P, b: P, c: P) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}
fn dist2(a: P, b: P) -> f64 {
    (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)
}
fn area(p: &[P]) -> f64 {
    (0..p.len())
        .map(|i| {
            let a = p[i];
            let b = p[(i + 1) % p.len()];
            a.0 * b.1 - a.1 * b.0
        })
        .sum()
}
fn on_segment(p: P, a: P, b: P) -> bool {
    p.0 >= a.0.min(b.0) - 1e-7
        && p.0 <= a.0.max(b.0) + 1e-7
        && p.1 >= a.1.min(b.1) - 1e-7
        && p.1 <= a.1.max(b.1) + 1e-7
        && cross(a, b, p).abs() < 1e-7
}
fn inside(p: P, poly: &[P]) -> bool {
    let mut inside = false;
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        if p.1 < a.1.min(b.1) - 1e-7 || p.1 > a.1.max(b.1) + 1e-7 {
            continue;
        }
        if on_segment(p, a, b) {
            return true;
        }
        if (a.1 > p.1) != (b.1 > p.1) && p.0 < a.0 + (p.1 - a.1) * (b.0 - a.0) / (b.1 - a.1) {
            inside = !inside;
        }
    }
    inside
}
fn proper_cross(a: P, b: P, c: P, d: P) -> bool {
    // Reject disjoint bounds before normalized orientation tests.
    if a.0.max(b.0) < c.0.min(d.0)
        || c.0.max(d.0) < a.0.min(b.0)
        || a.1.max(b.1) < c.1.min(d.1)
        || c.1.max(d.1) < a.1.min(b.1)
    {
        return false;
    }
    let ab = dist2(a, b).sqrt().max(1e-20);
    let cd = dist2(c, d).sqrt().max(1e-20);
    let opposite = |x: f64, y: f64| (x < -1e-7 && y > 1e-7) || (y < -1e-7 && x > 1e-7);
    opposite(cross(a, b, c) / ab, cross(a, b, d) / ab)
        && opposite(cross(c, d, a) / cd, cross(c, d, b) / cd)
}
fn contour_conflict(a: P, b: P, c: P, d: P) -> bool {
    if proper_cross(a, b, c, d) {
        return true;
    }
    // A T junction or overlapping collinear span can overlap two silhouettes
    // without a strict crossing. Even/odd face classification would cut a hole.
    let interior =
        |p: P, u: P, v: P| on_segment(p, u, v) && dist2(p, u) > 1e-12 && dist2(p, v) > 1e-12;
    interior(a, c, d)
        || interior(b, c, d)
        || interior(c, a, b)
        || interior(d, a, b)
        || (dist2(a, c) < 1e-12 && dist2(b, d) < 1e-12)
        || (dist2(a, d) < 1e-12 && dist2(b, c) < 1e-12)
}
fn valid(raw: &[P], candidate: &[P]) -> bool {
    if candidate.len() < 3 || area(raw) * area(candidate) <= 0.0 {
        return false;
    }
    for i in 0..candidate.len() {
        let a = candidate[i];
        let b = candidate[(i + 1) % candidate.len()];
        for j in i + 1..candidate.len() {
            if proper_cross(a, b, candidate[j], candidate[(j + 1) % candidate.len()]) {
                return false;
            }
        }
        for j in 0..raw.len() {
            if proper_cross(a, b, raw[j], raw[(j + 1) % raw.len()]) {
                return false;
            }
        }
    }
    let (inner, outer) = if area(raw) > 0.0 {
        (raw, candidate)
    } else {
        (candidate, raw)
    };
    inner.iter().enumerate().all(|(i, &p)| {
        let q = inner[(i + 1) % inner.len()];
        inside(p, outer) && inside(((p.0 + q.0) * 0.5, (p.1 + q.1) * 0.5), outer)
    })
}
fn approximate(raw: &[P], epsilon: f64, blockers: &[P]) -> Vec<P> {
    let n = raw.len();
    let split = (1..n)
        .max_by(|&a, &b| dist2(raw[0], raw[a]).total_cmp(&dist2(raw[0], raw[b])))
        .unwrap();
    let mut closed = raw.to_vec();
    closed.push(raw[0]);
    let mut keep = vec![false; n + 1];
    keep[0] = true;
    keep[split] = true;
    let mut stack = vec![(0, split), (split, n)];
    while let Some((a, b)) = stack.pop() {
        if b <= a + 1 {
            continue;
        }
        let u = closed[a];
        let v = closed[b];
        let length = dist2(u, v);
        // Allow more simplification on curved spans at low density. Retain a
        // 15% chord cap to keep enough curvature for stable pin deformation;
        // supporting lines and area/topology checks protect opaque coverage.
        let mut worst = epsilon.min((length.sqrt() * 0.15).max(1.25)).powi(2);
        let mut index = None;
        for (i, &p) in closed.iter().enumerate().take(b).skip(a + 1) {
            let t = if length == 0.0 {
                0.0
            } else {
                (((p.0 - u.0) * (v.0 - u.0) + (p.1 - u.1) * (v.1 - u.1)) / length).clamp(0.0, 1.0)
            };
            let d = dist2(p, (u.0 + t * (v.0 - u.0), u.1 + t * (v.1 - u.1)));
            if d > worst {
                worst = d;
                index = Some(i);
            }
        }
        if let Some(i) = index {
            keep[i] = true;
            stack.push((a, i));
            stack.push((i, b));
        }
    }
    let outer = area(raw) > 0.0;
    for _ in 0..32 {
        let selected = (0..n).filter(|&i| keep[i]).collect::<Vec<_>>();
        let candidate = supporting_polygon(&closed, &selected, epsilon);
        let blocked = blockers
            .iter()
            .copied()
            .filter(|&p| inside(p, &candidate) && !inside(p, raw))
            .collect::<Vec<_>>();
        if valid(raw, &candidate) && blocked.is_empty() {
            return candidate;
        }
        // Repair only the spans touching an invalid corner. A finger notch or
        // a single raster stair must not reduce tolerance around the whole body.
        let mut changed = false;
        // Preserve the local gap to another contour without lowering the error
        // tolerance around the entire silhouette (a one-pixel island can cause
        // thousands of unrelated raster vertices otherwise).
        for p in blocked {
            let j = (0..n)
                .min_by(|&a, &b| dist2(raw[a], p).total_cmp(&dist2(raw[b], p)))
                .unwrap();
            for k in [
                (j + n - 2) % n,
                (j + n - 1) % n,
                j,
                (j + 1) % n,
                (j + 2) % n,
            ] {
                if !keep[k] {
                    keep[k] = true;
                    changed = true;
                }
            }
        }
        for i in 0..candidate.len() {
            for j in i + 1..candidate.len() {
                if proper_cross(
                    candidate[i],
                    candidate[(i + 1) % candidate.len()],
                    candidate[j],
                    candidate[(j + 1) % candidate.len()],
                ) {
                    for p in [
                        candidate[i],
                        candidate[(i + 1) % candidate.len()],
                        candidate[j],
                        candidate[(j + 1) % candidate.len()],
                    ] {
                        let nearest = (0..selected.len())
                            .min_by(|&a, &b| {
                                dist2(raw[selected[a]], p).total_cmp(&dist2(raw[selected[b]], p))
                            })
                            .unwrap();
                        for k in [(nearest + selected.len() - 1) % selected.len(), nearest] {
                            let a = selected[k];
                            let b = if k + 1 < selected.len() {
                                selected[k + 1]
                            } else {
                                n
                            };
                            if b > a + 1 {
                                keep[(a + b) / 2] = true;
                                changed = true;
                            }
                        }
                    }
                }
            }
        }
        for j in 0..n {
            let p = raw[j];
            let q = raw[(j + 1) % n];
            let crossed = candidate
                .iter()
                .enumerate()
                .any(|(i, &a)| proper_cross(a, candidate[(i + 1) % candidate.len()], p, q));
            let outside = outer
                && (!inside(p, &candidate)
                    || !inside(((p.0 + q.0) * 0.5, (p.1 + q.1) * 0.5), &candidate));
            if crossed || outside {
                for k in [j, (j + 1) % n] {
                    if !keep[k] {
                        keep[k] = true;
                        changed = true;
                    }
                }
            }
        }
        if !changed {
            return candidate;
        }
    }
    supporting_polygon(
        &closed,
        &(0..n).filter(|&i| keep[i]).collect::<Vec<_>>(),
        epsilon,
    )
}
fn supporting_polygon(closed: &[P], selected: &[usize], epsilon: f64) -> Vec<P> {
    let n = closed.len() - 1;
    if selected.len() < 3 {
        return closed[..n].to_vec();
    }
    let mut lines = Vec::new();
    for k in 0..selected.len() {
        let a = selected[k];
        let b = if k + 1 < selected.len() {
            selected[k + 1]
        } else {
            n
        };
        let p = closed[a];
        let q = closed[b];
        let length = dist2(p, q).sqrt();
        let normal = ((q.1 - p.1) / length, (p.0 - q.0) / length);
        let offset = closed[a..=b]
            .iter()
            .map(|r| (r.0 - p.0) * normal.0 + (r.1 - p.1) * normal.1)
            .fold(0.0_f64, f64::max);
        lines.push((
            (p.0 + normal.0 * offset, p.1 + normal.1 * offset),
            (q.0 + normal.0 * offset, q.1 + normal.1 * offset),
        ));
    }
    let mut out = Vec::new();
    for i in 0..lines.len() {
        let (a, b) = lines[(i + lines.len() - 1) % lines.len()];
        let (c, d) = lines[i];
        let u = (b.0 - a.0, b.1 - a.1);
        let v = (d.0 - c.0, d.1 - c.1);
        let denominator = u.0 * v.1 - u.1 * v.0;
        if denominator.abs() > 1e-10 {
            let t = ((c.0 - a.0) * v.1 - (c.1 - a.1) * v.0) / denominator;
            let p = (a.0 + t * u.0, a.1 + t * u.1);
            if dist2(p, closed[selected[i]]) <= (epsilon * 2.0 + 1e-5).powi(2) {
                out.push(p);
                continue;
            }
        }
        out.push(b);
        if dist2(b, c) > 1e-14 {
            out.push(c);
        }
    }
    out.dedup_by(|a, b| dist2(*a, *b) < 1e-14);
    out
}
pub fn simplify(
    points: &[(f32, f32)],
    epsilon: f32,
    spacing: f32,
    area_budget: f64,
) -> Vec<(f32, f32)> {
    simplify_avoiding(points, epsilon, spacing, area_budget, &[])
}
fn simplify_avoiding(
    points: &[(f32, f32)],
    epsilon: f32,
    spacing: f32,
    area_budget: f64,
    blockers: &[P],
) -> Vec<(f32, f32)> {
    if points.len() < 3 {
        return points.to_vec();
    }
    let raw = points
        .iter()
        .enumerate()
        // Straight raster runs need only their endpoints. This is exact, and
        // prevents quadratic topology checks on thousands of redundant pixels.
        .filter(|&(i, &p)| {
            let a = points[(i + points.len() - 1) % points.len()];
            let b = points[(i + 1) % points.len()];
            let (a, p, b) = (
                (a.0 as f64, a.1 as f64),
                (p.0 as f64, p.1 as f64),
                (b.0 as f64, b.1 as f64),
            );
            cross(a, p, b) != 0.0 || (p.0 - a.0) * (b.0 - p.0) + (p.1 - a.1) * (b.1 - p.1) <= 0.0
        })
        .map(|(_, &(x, y))| (x as f64, y as f64))
        .collect::<Vec<_>>();
    let mut tolerance = if epsilon >= 4.0 {
        2.0_f64.powf((epsilon as f64).log2().floor())
    } else {
        (epsilon as f64 * 4.0).floor() * 0.25
    };
    let mut selected = raw.clone();
    // One pixel of perimeter area accounts for conservative stair quantization.
    // Without this allowance, a subpixel area budget restores every stair step.
    let raster_allowance = (0..raw.len())
        .map(|i| dist2(raw[i], raw[(i + 1) % raw.len()]).sqrt())
        .sum::<f64>()
        * 2.0;
    let mut unique = std::collections::HashSet::new();
    let touching = raw
        .iter()
        .any(|p| !unique.insert((p.0.to_bits(), p.1.to_bits())));
    let contracted = if touching {
        raw.clone()
    } else {
        contract(
            &raw,
            epsilon as f64,
            (area(&raw).abs() * area_budget).max(raster_allowance),
            blockers,
        )
    };
    // Repair this contour only; never globally revert the entire mesh to pixels.
    for _ in 0..24 {
        if touching {
            break;
        }
        let candidate = approximate(&raw, tolerance, blockers);
        if (area(&candidate) - area(&raw)).abs()
            <= (area(&raw).abs() * area_budget).max(raster_allowance) + 1e-8
            && valid(&raw, &candidate)
            && !blockers
                .iter()
                .any(|&p| inside(p, &candidate) && !inside(p, &raw))
        {
            selected = candidate;
            break;
        }
        // A shared tolerance ladder prevents density-dependent starting values
        // from skipping a better approximation under the same area allowance.
        tolerance = if tolerance > 4.0 {
            tolerance * 0.5
        } else {
            (tolerance - 0.25).max(0.0)
        };
    }
    if contracted.len() < selected.len() && valid(&raw, &contracted) {
        selected = contracted;
    }
    // Repair may restore short raster steps beside an otherwise simplified
    // span. Contract those steps on the resulting polygon as well, within the
    // original area allowance, rather than triangulating each step as a sliver.
    let budget = (area(&raw).abs() * area_budget).max(raster_allowance);
    let remaining = (budget - (area(&selected) - area(&raw)).abs()).max(0.0);
    let cleaned = if touching {
        selected.clone()
    } else {
        contract(&selected, epsilon as f64, remaining, blockers)
    };
    if cleaned.len() < selected.len() && valid(&raw, &cleaned) {
        selected = cleaned;
    }
    // Reserve a small local allowance for sub-spacing contour edges. Spending
    // it only on short edges preserves broad articulation seams such as necks.
    if !touching && area_budget >= 0.025 {
        let cleaned = contract_short_edges(
            &selected,
            epsilon as f64,
            area(&raw).abs() * area_budget * 0.2,
            blockers,
            epsilon as f64 * 0.75,
        );
        if cleaned.len() < selected.len() && valid(&raw, &cleaned) {
            selected = cleaned;
        }
    }
    let mut out = Vec::new();
    for i in 0..selected.len() {
        let a = selected[i];
        let b = selected[(i + 1) % selected.len()];
        let count = (dist2(a, b).sqrt() / spacing.max(1.0) as f64)
            .ceil()
            .max(1.0) as usize;
        for j in 0..count {
            let t = j as f64 / count as f64;
            out.push((
                (a.0 + t * (b.0 - a.0)) as f32,
                (a.1 + t * (b.1 - a.1)) as f32,
            ));
        }
    }
    out
}

pub fn contours(
    raw: &[Vec<(f32, f32)>],
    epsilon: f32,
    spacing: f32,
    area_budget: f64,
) -> Vec<Vec<(f32, f32)>> {
    let mut tolerances = vec![epsilon; raw.len()];
    let mut result = raw
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let blockers = raw
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .flat_map(|(_, c)| {
                    c.iter().enumerate().flat_map(|(k, &(x, y))| {
                        let q = c[(k + 1) % c.len()];
                        [
                            (x as f64, y as f64),
                            ((x as f64 + q.0 as f64) * 0.5, (y as f64 + q.1 as f64) * 0.5),
                        ]
                    })
                })
                .collect::<Vec<_>>();
            simplify_avoiding(p, epsilon, spacing, area_budget, &blockers)
        })
        .collect::<Vec<_>>();
    for attempt in 0..12 {
        let mut repair = vec![false; raw.len()];
        for a in 0..result.len() {
            for b in a + 1..result.len() {
                let overlaps = result[a].iter().enumerate().any(|(i, &p)| {
                    let q = result[a][(i + 1) % result[a].len()];
                    result[b].iter().enumerate().any(|(j, &r)| {
                        let s = result[b][(j + 1) % result[b].len()];
                        contour_conflict(
                            (p.0 as f64, p.1 as f64),
                            (q.0 as f64, q.1 as f64),
                            (r.0 as f64, r.1 as f64),
                            (s.0 as f64, s.1 as f64),
                        )
                    })
                });
                if overlaps {
                    repair[a] = true;
                    repair[b] = true;
                }
            }
        }
        if !repair.iter().any(|&r| r) {
            break;
        }
        for i in 0..raw.len() {
            if repair[i] {
                tolerances[i] = if attempt == 11 {
                    0.0
                } else {
                    tolerances[i] * 0.5
                };
                let blockers = raw
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != i)
                    .flat_map(|(_, c)| {
                        c.iter().enumerate().flat_map(|(k, &(x, y))| {
                            let q = c[(k + 1) % c.len()];
                            [
                                (x as f64, y as f64),
                                ((x as f64 + q.0 as f64) * 0.5, (y as f64 + q.1 as f64) * 0.5),
                            ]
                        })
                    })
                    .collect::<Vec<_>>();
                result[i] =
                    simplify_avoiding(&raw[i], tolerances[i], spacing, area_budget, &blockers);
            }
        }
    }
    for i in 0..result.len() {
        let polygon = result[i]
            .iter()
            .map(|&(x, y)| (x as f64, y as f64))
            .collect::<Vec<_>>();
        let blockers = result
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .flat_map(|(_, p)| p.iter().map(|&(x, y)| (x as f64, y as f64)))
            .collect::<Vec<_>>();
        let cleaned = clean_short_steps(
            &polygon,
            tolerances[i] as f64 * 0.75,
            tolerances[i] as f64,
            area(&polygon).abs() * area_budget * 0.2,
            &blockers,
        );
        let conflict = result
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .any(|(_, other)| {
                cleaned.iter().enumerate().any(|(a, &p)| {
                    other.iter().enumerate().any(|(b, &r)| {
                        let s = other[(b + 1) % other.len()];
                        contour_conflict(
                            p,
                            cleaned[(a + 1) % cleaned.len()],
                            (r.0 as f64, r.1 as f64),
                            (s.0 as f64, s.1 as f64),
                        )
                    })
                })
            });
        if !conflict {
            result[i] = cleaned
                .into_iter()
                .map(|(x, y)| (x as f32, y as f32))
                .collect();
        }
    }
    result
}
