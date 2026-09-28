use std::{
    cell::RefCell,
    cmp::Reverse,
    collections::{BTreeMap, BinaryHeap, HashMap},
    rc::Rc,
};

use aviutl2::anyhow::{self, Context as _};

const RELAXATION_STEPS: usize = 24;
const GRADIENT_CONTINUITY: f64 = 0.02;
const ROTATION_DIFFUSION: f64 = 0.001;
const STRAIN_WEIGHT: f64 = 64.;
const STRAIN_STEPS: usize = 24;
const POLAR_TRANSITION: f64 = 0.1;
use crate::linear::{PreparedSystem, SolveOptions};
#[cfg(test)]
thread_local! {
    static WORK: RefCell<[usize;3]> = const { RefCell::new([0;3]) };
    static SOLVE_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}
#[cfg(test)]
pub(crate) fn take_solve_calls() -> usize {
    SOLVE_CALLS.with(|n| n.replace(0))
}
#[cfg(test)]
pub(crate) fn take_work() -> [usize; 3] {
    let mut work = WORK.with(|w| std::mem::replace(&mut *w.borrow_mut(), [0; 3]));
    let linear = crate::linear::take_work();
    work[1] = linear[0];
    work[2] = linear[1];
    work
}

pub(crate) use crate::math::Point;

/// Planar ARAP initialized by a scale-free solve, with hard positional constraints.
/// A scale-free mesh solve selects a coherent rotation field before the rigid
/// reconstruction. Pin connectivity never supplies material rotations.
#[cfg(test)]
pub(crate) fn solve(
    vertices: &[f64],
    indices: &[i32],
    pin_source: &[f64],
    pin_destination: &[f64],
) -> anyhow::Result<Vec<Point>> {
    solve_with_poses(vertices, indices, pin_source, pin_destination, &[])
}

pub(crate) fn solve_with_poses(
    vertices: &[f64],
    indices: &[i32],
    pin_source: &[f64],
    pin_destination: &[f64],
    poses: &[crate::pose::Pose],
) -> anyhow::Result<Vec<Point>> {
    #[cfg(test)]
    SOLVE_CALLS.with(|n| n.set(n.get() + 1));
    let original = points(vertices)?;
    let triangles = triangles(indices, original.len())?;
    let sources = if pin_source.is_empty() {
        Vec::new()
    } else {
        points(pin_source)?
    };
    let destinations = if pin_destination.is_empty() {
        Vec::new()
    } else {
        points(pin_destination)?
    };
    anyhow::ensure!(
        sources.len() == destinations.len(),
        "pin arrays differ in length"
    );
    if sources == destinations && !poses.iter().any(|p| p.rotate || p.resize) {
        return Ok(original);
    }

    let adjacency = adjacency(&original, &triangles);
    // A single positional constraint has an exact translation solution. Keep
    // it exact (also in the legacy MLS entry point), without numerical drift.
    // Uncontrolled islands retain their own translation gauge.
    if sources.len() == 1 && poses.is_empty() {
        let translation = destinations[0].sub(sources[0]);
        let seed = nearest(&original, sources[0]);
        let mut result = original.clone();
        let mut visited = vec![false; original.len()];
        let mut stack = vec![seed];
        visited[seed] = true;
        while let Some(i) = stack.pop() {
            result[i] = original[i].add(translation);
            for &(j, _) in &adjacency[i] {
                if !visited[j] {
                    visited[j] = true;
                    stack.push(j);
                }
            }
        }
        return Ok(result);
    }
    let mut elements = elements(&original, &triangles);
    let widths = weight_thin_material(&original, &triangles, &mut elements);
    let thickness_weights = elements
        .iter()
        .zip(self::elements(&original, &triangles))
        .map(|(e, rest)| 4.0 * (3.0 - 6.0 * (rest.area / e.area).sqrt()).clamp(0.0, 1.0))
        .collect::<Vec<_>>();
    let scales = material_scales(&original, &adjacency, poses);
    let mut handles = handle_cells(&original, &triangles, &adjacency, poses);
    let patches = material_patches(
        &original, &triangles, &elements, &widths, &sources, poses, &scales,
    );
    let surface = element_adjacency(&elements, original.len());
    for patch in &patches {
        let ids = patch
            .samples
            .iter()
            .map(|s| s.0)
            .chain(std::iter::once(patch.center))
            .collect::<std::collections::BTreeSet<_>>();
        let mut mass = BTreeMap::<usize, f64>::new();
        for e in self::elements(&original, &triangles) {
            if e.vertices.iter().all(|i| ids.contains(i)) {
                for i in e.vertices {
                    *mass.entry(i).or_default() += e.area / 3.;
                }
            }
        }
        let total = mass.values().sum::<f64>();
        let centroid = mass.iter().fold(Point::default(), |p, (&i, &w)| {
            p.add(original[i].scale(w / total))
        });
        let radius2 = mass
            .iter()
            .map(|(&i, &w)| w * squared_distance(original[i], centroid))
            .sum::<f64>()
            / total;
        let frame = Rc::new(
            mass.iter()
                .map(|(&i, &w)| (i, w / total))
                .collect::<Vec<_>>(),
        );
        for (&i, &w) in &mass {
            handles.push(HandleCell {
                center: i,
                samples: mass
                    .range((i + 1)..)
                    .map(|(&j, &v)| (j, 8. * w * v / (total * radius2).max(1e-12)))
                    .collect(),
                angle: None,
                scale: patch.scale,
                frame: Some(frame.clone()),
            });
        }
    }
    let rigid_graph = combined_adjacency(&surface, &handles);
    let combined = gradient_regularization(
        &original,
        &triangles,
        &adjacency,
        &rigid_graph,
        GRADIENT_CONTINUITY,
    );
    let mut constraints = HashMap::<usize, Point>::new();
    for (&source, &destination) in sources.iter().zip(&destinations) {
        let vertex = nearest(&original, source);
        constraints.insert(vertex, original[vertex].add(destination.sub(source)));
    }
    anchor_unconstrained_components(&adjacency, &original, &mut constraints, poses);
    // A rigid region carrying a positional constraint must participate in the
    // solve so its exact target propagates across the material interface.
    // Free regions follow the solved field; projecting them on every step
    // would introduce artificial constraints and accumulate translational drift.
    let constrained_parts = patches
        .iter()
        .filter(|part| {
            constraints.contains_key(&part.center)
                || part.samples.iter().any(|s| constraints.contains_key(&s.0))
        })
        .cloned()
        .collect::<Vec<_>>();

    let mut deformed = original.clone();
    // Fit one rigid pose per connected component. This selects the natural
    // translation for a single handle and avoids the 180-degree identity trap.
    initialize_components(&adjacency, &original, &constraints, &mut deformed);
    for (&vertex, &target) in &constraints {
        deformed[vertex] = target;
    }

    let system = GlobalSystem::new(&combined, &constraints);
    let smoothing = SimilaritySmoother::new(&elements, &original);
    let mut first = conformal_pose(
        &original,
        &elements,
        &adjacency,
        &handles,
        &constraints,
        &deformed,
        &widths,
    );
    // Free explicit frames have no positional target to select a winding.
    // Normalizing the derivative and restoring lengths switches between two
    // distant minima when frames oppose each other. Keep the linear frame
    // reconstruction in this case: it is periodic and continuous in each
    // rotation, while equal frames still reproduce a rigid transformation.
    // Positional rigs and freely following starch handles keep the ARAP solve.
    if sources.is_empty() && !poses.is_empty() && poses.iter().all(|p| !p.moving && p.rotate) {
        return Ok(first);
    }
    // Both phases have the same fixed iteration budget for every target pose.
    // Only the first phase regularizes derivatives: the final reconstruction
    // permits a shoulder fold instead of flattening its material against a barrier.
    let protected = elements
        .iter()
        .map(|e| e.vertices.iter().any(|i| constraints.contains_key(i)))
        .collect::<Vec<_>>();
    let rigid_system = GlobalSystem::new(&rigid_graph, &constraints);
    let empty = vec![Vec::new(); original.len()];
    let unit = vec![(1., 0.); original.len()];
    for (current_system, strain_weight, count) in [
        (&system, 0., RELAXATION_STEPS),
        (&rigid_system, STRAIN_WEIGHT, STRAIN_STEPS),
    ] {
        let reference = material_similarities(&elements, &first, &smoothing);
        for step in 0..count {
            let strain_weight = strain_weight * (step + 1) as f64 / count as f64;
            let mut current = material_similarities(&elements, &first, &smoothing);
            if strain_weight > 0. {
                for (z, &initial) in current.iter_mut().zip(&reference) {
                    *z = z.add(initial).scale(0.5);
                }
            }
            let rotations = &current;
            let mut rhs = element_rhs(
                &elements,
                &first,
                rotations,
                &scales,
                strain_weight,
                &protected,
                &thickness_weights,
            );
            for p in &mut rhs {
                *p = p.scale(1. / (1. + strain_weight));
            }
            for h in &handles {
                let rotation = h.rotation(&original, &first);
                for &(j, w) in &h.samples {
                    let edge =
                        rotate(original[j].sub(original[h.center]), rotation).scale(h.scale * w);
                    rhs[j] = rhs[j].add(edge);
                    rhs[h.center] = rhs[h.center].sub(edge);
                }
            }
            let next =
                current_system.step_extra(&original, &empty, &unit, &constraints, &first, &rhs);
            if strain_weight > 0.0 {
                for (p, q) in first.iter_mut().zip(next) {
                    *p = p.scale(0.75).add(q.scale(0.25));
                }
            } else {
                first = next;
            }
            if strain_weight > 0.0 {
                project_material_parts(
                    &original,
                    &triangles,
                    &constraints,
                    &constrained_parts,
                    &mut first,
                );
            }
        }
    }
    project_material_parts(&original, &triangles, &constraints, &patches, &mut first);
    Ok(first)
}

fn cross_section(boundary: &[(Point, Point)], center: Point) -> (f64, Point) {
    let mut transverse = Point { x: 1., y: 0. };
    let cross = |a: Point, b: Point| a.x * b.y - a.y * b.x;
    let mut width = f64::INFINITY;
    for k in 0..16 {
        let angle = k as f64 * std::f64::consts::PI / 16.;
        let v = Point {
            x: angle.cos(),
            y: angle.sin(),
        };
        let (mut positive, mut negative) = (f64::INFINITY, f64::INFINITY);
        for &(a, b) in boundary {
            let edge = b.sub(a);
            let denom = cross(v, edge);
            if denom.abs() < 1e-12 {
                continue;
            }
            let t = cross(a.sub(center), edge) / denom;
            let u = cross(a.sub(center), v) / denom;
            if !(0. ..=1.).contains(&u) {
                continue;
            }
            if t >= 0. {
                positive = positive.min(t);
            } else {
                negative = negative.min(-t);
            }
        }
        if positive + negative < width {
            width = positive + negative;
            transverse = v;
        }
    }
    (width, transverse)
}

fn edge_multiplicities(
    triangles: impl IntoIterator<Item = [usize; 3]>,
) -> BTreeMap<(usize, usize), usize> {
    let mut counts = BTreeMap::new();
    for [a, b, c] in triangles {
        for (a, b) in [(a, b), (b, c), (c, a)] {
            *counts.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    counts
}

struct MaterialCache {
    original: Vec<Point>,
    triangles: Vec<[usize; 3]>,
    areas: Vec<f64>,
    normals: Vec<Point>,
    widths: Vec<f64>,
}

thread_local! {static MATERIAL:RefCell<Option<MaterialCache>>=const{RefCell::new(None)};}
// Equalize narrow and broad material stiffness using rest-domain cross sections.
// This does not change a pin's range and has no dependence on its destination.
fn weight_thin_material(
    original: &[Point],
    triangles: &[[usize; 3]],
    elements: &mut [Element],
) -> Vec<f64> {
    let hit = MATERIAL.with(|cache| {
        let cache = cache.borrow();
        if let Some(c) = cache.as_ref() {
            if c.original == original && c.triangles == triangles {
                for ((e, &a), &normal) in elements.iter_mut().zip(&c.areas).zip(&c.normals) {
                    e.area = a;
                    e.transverse = normal;
                }
                return Some(c.widths.clone());
            }
        }
        None
    });
    if let Some(widths) = hit {
        return widths;
    }
    let widths = compute_material_weights(original, triangles, elements);
    MATERIAL.with(|cache| {
        *cache.borrow_mut() = Some(MaterialCache {
            original: original.to_vec(),
            triangles: triangles.to_vec(),
            areas: elements.iter().map(|e| e.area).collect(),
            normals: elements.iter().map(|e| e.transverse).collect(),
            widths: widths.clone(),
        })
    });
    widths
}
fn compute_material_weights(
    original: &[Point],
    triangles: &[[usize; 3]],
    elements: &mut [Element],
) -> Vec<f64> {
    let edges = edge_multiplicities(triangles.iter().copied());
    let boundary = edges
        .into_iter()
        .filter(|(_, n)| *n == 1)
        .map(|((a, b), _)| (original[a], original[b]))
        .collect::<Vec<_>>();
    let widths = elements
        .iter_mut()
        .map(|e| {
            let center = e
                .vertices
                .iter()
                .fold(Point::default(), |q, &i| q.add(original[i]))
                .scale(1. / 3.);
            let (width, transverse) = cross_section(&boundary, center);
            e.transverse = transverse;
            width
        })
        .collect::<Vec<_>>();
    // Normalize within each connected material island.
    let mut component = (0..original.len()).collect::<Vec<_>>();
    fn root(ids: &mut [usize], i: usize) -> usize {
        if ids[i] != i {
            ids[i] = root(ids, ids[i]);
        }
        ids[i]
    }
    for &[a, b, c] in triangles {
        let a = root(&mut component, a);
        for i in [b, c] {
            let b = root(&mut component, i);
            component[b] = a;
        }
    }
    let mut maxima = vec![0.; original.len()];
    for (e, &width) in elements.iter().zip(&widths) {
        let id = root(&mut component, e.vertices[0]);
        if width.is_finite() {
            maxima[id] = f64::max(maxima[id], width);
        }
    }
    let mut vertex_widths = vec![0.; original.len()];
    for (e, &width) in elements.iter().zip(&widths) {
        for &i in &e.vertices {
            if width.is_finite() {
                vertex_widths[i] = f64::max(vertex_widths[i], width);
            }
        }
    }
    for (e, width) in elements.iter_mut().zip(widths) {
        let maximum = maxima[root(&mut component, e.vertices[0])];
        e.area *= if width.is_finite() && width > 1e-12 {
            (maximum / width).powi(2).clamp(1., 16.)
        } else {
            1.
        };
    }
    vertex_widths
}

fn project_material_parts(
    original: &[Point],
    triangles: &[[usize; 3]],
    constraints: &HashMap<usize, Point>,
    parts: &[HandleCell],
    output: &mut [Point],
) {
    let mut sums = vec![Point::default(); original.len()];
    let mut weights = vec![0.; original.len()];
    let edge_counts = edge_multiplicities(triangles.iter().copied());
    for part in parts {
        let mut vertices = part
            .samples
            .iter()
            .map(|s| s.0)
            .collect::<std::collections::BTreeSet<_>>();
        vertices.insert(part.center);
        let anchor = vertices
            .iter()
            .find(|i| constraints.contains_key(i))
            .copied();
        // Project the deformed material onto its closest rigid frame. Integrate
        // over triangle area so contour sampling and pin type cannot bias it.
        // A positional constraint fixes translation; it never prescribes angle.
        let faces = triangles
            .iter()
            .filter(|t| t.iter().all(|i| vertices.contains(i)))
            .copied()
            .collect::<Vec<_>>();
        let part_edges = edge_multiplicities(faces.iter().copied());
        let mut seam_mass = 0.;
        let mut seam_source = Point::default();
        let mut seam_target = Point::default();
        for (&(a, b), &count) in &part_edges {
            if count == 1 && edge_counts[&(a, b)] == 2 {
                let length = squared_distance(original[a], original[b]).sqrt();
                seam_mass += length;
                seam_source = seam_source.add(original[a].add(original[b]).scale(length * 0.5));
                seam_target = seam_target.add(output[a].add(output[b]).scale(length * 0.5));
            }
        }
        let seam = (seam_mass > 1e-12).then(|| {
            (
                seam_source.scale(1. / seam_mass),
                seam_target.scale(1. / seam_mass),
            )
        });
        // A nearby real control defines the joint pivot. Do not attach a free
        // part to a distant wrist merely because it is the nearest control.
        let joint = anchor
            .is_none()
            .then(|| {
                seam.and_then(|(p, _)| {
                    constraints
                        .keys()
                        .copied()
                        .filter(|&i| squared_distance(original[i], p) <= seam_mass * seam_mass)
                        .min_by(|&a, &b| {
                            squared_distance(original[a], p)
                                .total_cmp(&squared_distance(original[b], p))
                                .then(a.cmp(&b))
                        })
                })
            })
            .flatten();
        let mut samples = Vec::new();
        for [a, b, c] in faces {
            let u = original[b].sub(original[a]);
            let v = original[c].sub(original[a]);
            let area = (u.x * v.y - u.y * v.x).abs() * 0.5;
            for bary in [
                [2. / 3., 1. / 6., 1. / 6.],
                [1. / 6., 2. / 3., 1. / 6.],
                [1. / 6., 1. / 6., 2. / 3.],
            ] {
                let sample = |field: &[Point]| {
                    field[a]
                        .scale(bary[0])
                        .add(field[b].scale(bary[1]))
                        .add(field[c].scale(bary[2]))
                };
                samples.push((sample(original), sample(output), area / 3.));
            }
        }
        let mass = samples.iter().map(|s| s.2).sum::<f64>();
        if mass <= 1e-12 {
            continue;
        }
        let (p, q) = if let Some(i) = anchor.or(joint) {
            (original[i], constraints[&i])
        } else if let Some(seam) = seam {
            seam
        } else {
            (
                samples
                    .iter()
                    .fold(Point::default(), |p, s| p.add(s.0.scale(s.2 / mass))),
                samples
                    .iter()
                    .fold(Point::default(), |p, s| p.add(s.1.scale(s.2 / mass))),
            )
        };
        let mut z = Point::default();
        for (a, b, w) in samples {
            let a = a.sub(p);
            let b = b.sub(q);
            z.x += w * (a.x * b.x + a.y * b.y);
            z.y += w * (a.x * b.y - a.y * b.x);
        }
        if let Some(anchor) = joint {
            // Fit an unpinned part about its joint, not its displaced centroid.
            // The same frame controls both orientation and attachment position.
            // Inverse squared distance gives each control equal angular weight
            // and preserves rigid rotations without using bone connectivity.
            let mut control_rotation = Point::default();
            let mut control_count = 0;
            let mut controls = constraints.iter().collect::<Vec<_>>();
            controls.sort_by_key(|(i, _)| **i);
            for (&i, &target) in controls {
                if i == anchor {
                    continue;
                }
                let a = original[i].sub(original[anchor]);
                let b = target.sub(constraints[&anchor]);
                let length2 = a.x * a.x + a.y * a.y;
                if length2 > 1e-12 {
                    control_rotation.x += (a.x * b.x + a.y * b.y) / length2;
                    control_rotation.y += (a.x * b.y - a.y * b.x) / length2;
                    control_count += 1;
                }
            }
            if control_count > 0 && control_rotation.x.hypot(control_rotation.y) > 1e-12 {
                z = control_rotation;
            }
        }
        let norm = z.x.hypot(z.y);
        let rotation = if norm > 1e-12 {
            (z.x / norm, z.y / norm)
        } else {
            part.rotation(original, output)
        };
        for i in vertices {
            let target = q.add(rotate(original[i].sub(p), rotation).scale(part.scale));
            sums[i] = sums[i].add(target);
            weights[i] += 1.;
        }
    }
    for i in 0..output.len() {
        if weights[i] > 0. && !constraints.contains_key(&i) {
            output[i] = sums[i].scale(1. / weights[i]);
        }
    }
}

/// Rest-domain shape constraints used by the shared local/global solver.
fn material_patches(
    original: &[Point],
    triangles: &[[usize; 3]],
    weighted: &[Element],
    widths: &[f64],
    sources: &[Point],
    poses: &[crate::pose::Pose],
    scales: &[f64],
) -> Vec<HandleCell> {
    let raw = elements(original, triangles);
    let n = raw.len();
    let width = raw
        .iter()
        .zip(weighted)
        .map(|(a, b)| (a.area / b.area).sqrt())
        .collect::<Vec<_>>();
    let counts = edge_multiplicities(raw.iter().map(|e| e.vertices));
    let boundary = counts
        .iter()
        .filter(|(_, n)| **n == 1)
        .map(|(&(a, b), _)| (original[a], original[b]))
        .collect::<Vec<_>>();
    let maximum = widths.iter().copied().fold(0., f64::max);
    let thickness = |p: Point| cross_section(&boundary, p).0 / maximum;
    let mut incident = BTreeMap::new();
    let mut edges = Vec::new();
    for (i, e) in raw.iter().enumerate() {
        let [a, b, c] = e.vertices;
        for (a, b) in [(a, b), (b, c), (c, a)] {
            let key = (a.min(b), a.max(b));
            if let Some(&j) = incident.get(&key) {
                let saddle = thickness(original[a].add(original[b]).scale(0.5));
                edges.push((width[i].min(width[j]).min(saddle), i, j));
            } else {
                incident.insert(key, i);
            }
        }
    }
    edges.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    let mut ids = (0..n).collect::<Vec<_>>();
    let mut peaks = width.clone();
    fn root(ids: &mut [usize], i: usize) -> usize {
        if ids[i] != i {
            ids[i] = root(ids, ids[i]);
        }
        ids[i]
    }
    let links = edges.iter().map(|e| (e.1, e.2)).collect::<Vec<_>>();
    for (saddle, a, b) in edges {
        let a = root(&mut ids, a);
        let b = root(&mut ids, b);
        if a != b && peaks[a].min(peaks[b]) <= 1.25 * saddle {
            ids[b] = a;
            peaks[a] = peaks[a].max(peaks[b]);
        }
    }
    let mut groups = BTreeMap::<usize, Vec<usize>>::new();
    for i in 0..n {
        let id = root(&mut ids, i);
        groups.entry(id).or_default().push(i);
    }
    // Peripheral material too thin to become a rigid part still belongs to its
    // sole substantial neighbor. Use material width as well as face count:
    // an ear with three triangles must not stop following the head merely
    // because a different triangulation split its two-triangle fringe.
    let mut fringe_links = BTreeMap::<usize, Vec<usize>>::new();
    for &(a, b) in &links {
        let (a, b) = (root(&mut ids, a), root(&mut ids, b));
        if a != b {
            fringe_links.entry(a).or_default().push(b);
            fringe_links.entry(b).or_default().push(a);
        }
    }
    let mut visited = std::collections::BTreeSet::new();
    let is_fringe = |id: usize| {
        let area = groups[&id].iter().map(|&i| raw[i].area).sum::<f64>();
        // A large triangle can give a tiny peripheral lobe a high peak width.
        // It needs enough material area to support that width as a real part.
        groups[&id].len() < 3 || peaks[id] < 0.5 || area < (peaks[id] * maximum).powi(2) * 0.15
    };
    for (&id, _) in &groups {
        if !is_fringe(id) || !visited.insert(id) {
            continue;
        }
        let mut fringe = vec![id];
        let mut neighbors = std::collections::BTreeSet::new();
        let mut cursor = 0;
        while cursor < fringe.len() {
            let i = fringe[cursor];
            cursor += 1;
            for &j in fringe_links.get(&i).into_iter().flatten() {
                if !is_fringe(j) {
                    neighbors.insert(j);
                } else if visited.insert(j) {
                    fringe.push(j);
                }
            }
        }
        // A fringe bridging two parts remains flexible; never merge a joint.
        if neighbors.len() == 1 {
            let neighbor = *neighbors.first().unwrap();
            for i in fringe {
                ids[i] = neighbor;
            }
        }
    }
    groups.clear();
    for i in 0..n {
        let id = root(&mut ids, i);
        groups.entry(id).or_default().push(i);
    }
    let mut controls = BTreeMap::new();
    for &p in sources {
        controls.insert(nearest(original, p), false);
    }
    for p in poses.iter().filter(|p| p.moving || p.rotate || p.resize) {
        let i = nearest(
            original,
            Point {
                x: p.source.0,
                y: p.source.1,
            },
        );
        *controls.entry(i).or_default() |= p.rotate || (p.scale - 1.).abs() > 1e-8;
    }
    let mut owners = BTreeMap::<usize, Vec<(usize, bool)>>::new();
    for (vertex, explicit) in controls {
        let mut mass = BTreeMap::<usize, f64>::new();
        for (i, e) in raw
            .iter()
            .enumerate()
            .filter(|(_, e)| e.vertices.contains(&vertex))
        {
            *mass.entry(root(&mut ids, i)).or_default() += e.area;
        }
        if let Some((&id, _)) = mass.iter().max_by(|a, b| a.1.total_cmp(b.1)) {
            owners.entry(id).or_default().push((vertex, explicit));
        }
    }
    let mut patches = Vec::new();
    for (id, faces) in groups {
        let controls = owners.get(&id).map(Vec::as_slice).unwrap_or(&[]);
        #[cfg(test)]
        if std::env::var_os("PUPPET_TRACE_PARTS").is_some() {
            let min_y = faces
                .iter()
                .flat_map(|&f| raw[f].vertices)
                .map(|i| original[i].y)
                .fold(f64::INFINITY, f64::min);
            let max_y = faces
                .iter()
                .flat_map(|&f| raw[f].vertices)
                .map(|i| original[i].y)
                .fold(f64::NEG_INFINITY, f64::max);
            eprintln!(
                "part {id} faces {} peak {} y {min_y}..{max_y} controls {:?}",
                faces.len(),
                peaks[id],
                controls
                    .iter()
                    .map(|c| (original[c.0], c.1))
                    .collect::<Vec<_>>()
            );
        }
        if faces.len() < 3 || peaks[id] < 0.5 || controls.len() > 1 || controls.iter().any(|c| c.1)
        {
            continue;
        }
        // Count connected attachment seams, not neighboring watershed labels:
        // one neck can touch both halves of a split torso region.
        let edge_counts = edge_multiplicities(faces.iter().map(|&i| raw[i].vertices));
        let mut seam = BTreeMap::<usize, Vec<usize>>::new();
        for ((a, b), n) in edge_counts {
            if n == 1 && counts[&(a, b)] > 1 {
                seam.entry(a).or_default().push(b);
                seam.entry(b).or_default().push(a);
            }
        }
        let mut visited = std::collections::BTreeSet::new();
        let mut seams = 0;
        for &i in seam.keys() {
            if !visited.insert(i) {
                continue;
            }
            seams += 1;
            let mut stack = vec![i];
            while let Some(j) = stack.pop() {
                for &k in &seam[&j] {
                    if visited.insert(k) {
                        stack.push(k);
                    }
                }
            }
        }
        #[cfg(test)]
        if std::env::var_os("PUPPET_TRACE_PARTS").is_some() {
            eprintln!("part {id} seams {seams}");
        }
        if seams > 1 {
            continue;
        }
        let mut masses = BTreeMap::<usize, f64>::new();
        for i in faces {
            for v in raw[i].vertices {
                *masses.entry(v).or_default() += raw[i].area / 3.;
            }
        }
        let mass = masses.values().sum::<f64>();
        let centroid = masses.iter().fold(Point::default(), |q, (&i, &w)| {
            q.add(original[i].scale(w / mass))
        });
        let center = *masses
            .keys()
            .min_by(|&&a, &&b| {
                squared_distance(original[a], centroid)
                    .total_cmp(&squared_distance(original[b], centroid))
            })
            .unwrap();
        let radius2 = masses
            .iter()
            .map(|(&i, &w)| w * squared_distance(original[i], original[center]))
            .sum::<f64>()
            / mass;
        if radius2 < 1e-12 {
            continue;
        }
        let samples = masses
            .into_iter()
            .filter(|(i, _)| *i != center)
            .map(|(i, w)| (i, 8. * w / radius2))
            .collect();
        patches.push(HandleCell {
            frame: None,
            center,
            samples,
            angle: None,
            scale: scales[center],
        });
    }
    patches
}

// User scale is a differential material property. Harmonic extension carries
// it along unconstrained extremities instead of expanding only a shoulder disk.
fn material_scales(
    original: &[Point],
    adjacency: &[Vec<(usize, f64)>],
    poses: &[crate::pose::Pose],
) -> Vec<f64> {
    if !poses.iter().any(|p| (p.scale - 1.0).abs() > 1e-12) {
        return vec![1.0; original.len()];
    }
    let mut fixed = HashMap::new();
    for p in poses.iter().filter(|p| p.moving || p.resize) {
        fixed.insert(
            nearest(
                original,
                Point {
                    x: p.source.0,
                    y: p.source.1,
                },
            ),
            Point {
                x: p.scale,
                y: p.scale,
            },
        );
    }
    let ones = vec![Point { x: 1.0, y: 1.0 }; original.len()];
    anchor_unconstrained_components(adjacency, &ones, &mut fixed, &[]);
    let zero = vec![Point::default(); original.len()];
    GlobalSystem::new(adjacency, &fixed)
        .step_extra(
            &zero,
            adjacency,
            &vec![(1.0, 0.0); original.len()],
            &fixed,
            &ones,
            &zero,
        )
        .iter()
        .map(|p| p.x.max(0.0))
        .collect()
}

// Additional overlapping ARAP cells, one per finite material handle:
// E_h = sum_j w_hj |(q_j-q_h) - s_h R_h (p_j-p_h)|^2.
// They constrain relative shape, never freeze a translated disk of vertices.
// R_h is minimized each local step unless the user explicitly controls it.
#[derive(Clone)]
struct HandleCell {
    frame: Option<Rc<Vec<(usize, f64)>>>,
    center: usize,
    samples: Vec<(usize, f64)>,
    angle: Option<f64>,
    scale: f64,
}
impl HandleCell {
    fn rotation(&self, original: &[Point], deformed: &[Point]) -> (f64, f64) {
        if let Some(angle) = self.angle {
            return (angle.cos(), angle.sin());
        }
        let (mut c, mut s) = (0.0, 0.0);
        if let Some(frame) = &self.frame {
            let p = frame
                .iter()
                .fold(Point::default(), |p, &(i, w)| p.add(original[i].scale(w)));
            let q = frame
                .iter()
                .fold(Point::default(), |p, &(i, w)| p.add(deformed[i].scale(w)));
            for &(i, w) in frame.iter() {
                let a = original[i].sub(p);
                let b = deformed[i].sub(q);
                c += w * (a.x * b.x + a.y * b.y);
                s += w * (a.x * b.y - a.y * b.x);
            }
            let norm = c.hypot(s);
            return if norm > 1e-12 {
                (c / norm, s / norm)
            } else {
                (1., 0.)
            };
        }
        for &(j, w) in &self.samples {
            let p = original[j].sub(original[self.center]);
            let q = deformed[j].sub(deformed[self.center]);
            c += w * (p.x * q.x + p.y * q.y);
            s += w * (p.x * q.y - p.y * q.x);
        }
        unit_rotation(c, s)
    }
}
fn handle_cells(
    original: &[Point],
    triangles: &[[usize; 3]],
    adjacency: &[Vec<(usize, f64)>],
    poses: &[crate::pose::Pose],
) -> Vec<HandleCell> {
    let mut areas = vec![0.0; original.len()];
    for &[a, b, c] in triangles {
        let u = original[b].sub(original[a]);
        let v = original[c].sub(original[a]);
        let area = (u.x * v.y - u.y * v.x).abs() / 6.0;
        for i in [a, b, c] {
            areas[i] += area;
        }
    }
    poses
        .iter()
        .filter(|p| p.moving || p.rotate || p.resize)
        .map(|p| {
            let center = nearest(
                original,
                Point {
                    x: p.source.0,
                    y: p.source.1,
                },
            );
            let mut distances = vec![f64::INFINITY; original.len()];
            distances[center] = 0.0;
            let mut queue = vec![center];
            let mut cursor = 0;
            while cursor < queue.len() {
                let i = queue[cursor];
                cursor += 1;
                for &(j, _) in &adjacency[i] {
                    let d = distances[i] + squared_distance(original[i], original[j]).sqrt();
                    if d < 2.0 * p.radius && d < distances[j] - 1e-10 {
                        distances[j] = d;
                        queue.push(j);
                    }
                }
            }
            let mut samples: Vec<_> = distances
                .iter()
                .enumerate()
                .filter_map(|(j, &d)| {
                    let mut w = crate::pose::weight(d, p.radius) * areas[j];
                    for other in poses.iter().filter(|o| {
                        o.moving && (o.source.0 - p.source.0).hypot(o.source.1 - p.source.1) > 1e-8
                    }) {
                        let gap =
                            (original[j].x - other.source.0).hypot(original[j].y - other.source.1);
                        let t = (gap / (d + 1e-12) - 1.0).clamp(0.0, 1.0);
                        w *= t * t * (3.0 - 2.0 * t);
                    }
                    (j != center && w > 0.0).then_some((j, w))
                })
                .collect();
            let coarse = samples.len() < 2;
            if coarse {
                // A coarse mesh may have no vertices inside the physical handle
                // radius. Use its existing incident material instead of adding
                // a tiny refinement ring or silently dropping the handle.
                samples = adjacency[center]
                    .iter()
                    .filter_map(|&(j, _)| {
                        let controlled = poses.iter().any(|other| {
                            other.moving
                                && (original[j].x - other.source.0)
                                    .hypot(original[j].y - other.source.1)
                                    < 1e-8
                        });
                        (!controlled).then_some((
                            j,
                            areas[j] / squared_distance(original[j], original[center]).max(1e-12),
                        ))
                    })
                    .collect();
            }
            let sum = samples.iter().map(|s| s.1).sum::<f64>().max(1e-20);
            for sample in &mut samples {
                sample.1 *= if p.moving && !p.rotate {
                    64.0 / sum
                } else if p.rotate {
                    512.0 / sum
                } else {
                    384.0 / sum
                };
            }
            HandleCell {
                frame: None,
                center,
                samples,
                angle: p.rotate.then_some(p.natural_angle + p.angle),
                scale: p.scale,
            }
        })
        .collect()
}
fn combined_adjacency(
    base: &[Vec<(usize, f64)>],
    handles: &[HandleCell],
) -> Vec<Vec<(usize, f64)>> {
    let mut edges: Vec<BTreeMap<usize, f64>> = base
        .iter()
        .map(|row| row.iter().copied().collect())
        .collect();
    for h in handles {
        for &(j, w) in &h.samples {
            *edges[h.center].entry(j).or_default() += w;
            *edges[j].entry(h.center).or_default() += w;
        }
    }
    edges
        .into_iter()
        .map(|row| row.into_iter().collect())
        .collect()
}

fn points(values: &[f64]) -> anyhow::Result<Vec<Point>> {
    anyhow::ensure!(
        values.len() >= 2 && values.len().is_multiple_of(2),
        "invalid point array"
    );
    anyhow::ensure!(
        values.iter().all(|v| v.is_finite()),
        "point array contains non-finite values"
    );
    Ok(values
        .chunks_exact(2)
        .map(|v| Point { x: v[0], y: v[1] })
        .collect())
}

fn triangles(values: &[i32], vertex_count: usize) -> anyhow::Result<Vec<[usize; 3]>> {
    anyhow::ensure!(
        values.len() >= 3 && values.len().is_multiple_of(3),
        "invalid index array"
    );
    values
        .chunks_exact(3)
        .map(|t| {
            let result = [
                usize::try_from(t[0]).context("negative index")?,
                usize::try_from(t[1]).context("negative index")?,
                usize::try_from(t[2]).context("negative index")?,
            ];
            anyhow::ensure!(
                result.iter().all(|&i| i < vertex_count),
                "index out of bounds"
            );
            Ok(result)
        })
        .collect()
}

fn adjacency(original: &[Point], triangles: &[[usize; 3]]) -> Vec<Vec<(usize, f64)>> {
    let mut edges = vec![BTreeMap::<usize, f64>::new(); original.len()];
    for &[a, b, c] in triangles {
        for (i, j, k) in [(a, b, c), (b, c, a), (c, a, b)] {
            let u = original[i].sub(original[k]);
            let v = original[j].sub(original[k]);
            let area = (u.x * v.y - u.y * v.x).abs();
            if area <= 1.0e-14 {
                continue;
            }
            // Accumulate both opposite angles before enforcing nonnegative
            // edge weights. Negative weights can give a boundary cell a
            // negative rigidity energy and rotate even an unchanged mesh.
            let weight = 0.5 * (u.x * v.x + u.y * v.y) / area;
            *edges[i].entry(j).or_default() += weight;
            *edges[j].entry(i).or_default() += weight;
        }
    }
    edges
        .into_iter()
        .map(|e| e.into_iter().map(|(j, w)| (j, w.max(0.0))).collect())
        .collect()
}

fn initialize_components(
    adjacency: &[Vec<(usize, f64)>],
    original: &[Point],
    constraints: &HashMap<usize, Point>,
    output: &mut [Point],
) {
    let mut visited = vec![false; original.len()];
    for start in 0..original.len() {
        if visited[start] {
            continue;
        }
        let mut component = vec![start];
        visited[start] = true;
        let mut cursor = 0;
        while cursor < component.len() {
            for &(j, _) in &adjacency[component[cursor]] {
                if !visited[j] {
                    visited[j] = true;
                    component.push(j);
                }
            }
            cursor += 1;
        }
        let handles = component
            .iter()
            .copied()
            .filter(|i| constraints.contains_key(i))
            .collect::<Vec<_>>();
        if handles.is_empty() {
            continue;
        }
        let mut p = Point::default();
        let mut q = Point::default();
        for &i in &handles {
            p = p.add(original[i]);
            q = q.add(constraints[&i]);
        }
        p = p.scale(1.0 / handles.len() as f64);
        q = q.scale(1.0 / handles.len() as f64);
        let (mut real, mut imaginary) = (0.0, 0.0);
        for &i in &handles {
            let a = original[i].sub(p);
            let b = constraints[&i].sub(q);
            real += a.x * b.x + a.y * b.y;
            imaginary += a.x * b.y - a.y * b.x;
        }
        let rotation = unit_rotation(real, imaginary);
        for i in component {
            output[i] = rotate(original[i].sub(p), rotation).add(q);
        }
    }
}

fn nearest(points: &[Point], target: Point) -> usize {
    points
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            squared_distance(**a, target).total_cmp(&squared_distance(**b, target))
        })
        .map(|(index, _)| index)
        .unwrap()
}

fn squared_distance(a: Point, b: Point) -> f64 {
    let d = a.sub(b);
    d.x * d.x + d.y * d.y
}

fn anchor_unconstrained_components(
    adjacency: &[Vec<(usize, f64)>],
    original: &[Point],
    constraints: &mut HashMap<usize, Point>,
    poses: &[crate::pose::Pose],
) {
    let mut visited = vec![false; original.len()];
    for start in 0..original.len() {
        if visited[start] {
            continue;
        }
        let mut stack = vec![start];
        let mut component = Vec::new();
        let mut constrained = false;
        visited[start] = true;
        while let Some(vertex) = stack.pop() {
            component.push(vertex);
            constrained |= constraints.contains_key(&vertex);
            for &(next, _) in &adjacency[vertex] {
                if !visited[next] {
                    visited[next] = true;
                    stack.push(next);
                }
            }
        }
        if !constrained {
            let frame = poses.iter().filter(|p| p.rotate || p.resize).find_map(|p| {
                let i = nearest(
                    original,
                    Point {
                        x: p.source.0,
                        y: p.source.1,
                    },
                );
                component.contains(&i).then_some((i, p))
            });
            if let Some((i, p)) = frame {
                constraints.insert(
                    i,
                    original[i].add(Point {
                        x: p.center.0 - p.source.0,
                        y: p.center.1 - p.source.1,
                    }),
                );
            } else {
                constraints.insert(component[0], original[component[0]]);
            }
        }
    }
}

#[cfg(test)]
fn local_step(
    original: &[Point],
    deformed: &[Point],
    adjacency: &[Vec<(usize, f64)>],
) -> Vec<(f64, f64)> {
    (0..original.len())
        .map(|i| {
            let mut cosine = 0.0;
            let mut sine = 0.0;
            for &(j, weight) in &adjacency[i] {
                let p = original[i].sub(original[j]);
                let q = deformed[i].sub(deformed[j]);
                cosine += weight * (p.x * q.x + p.y * q.y);
                sine += weight * (p.x * q.y - p.y * q.x);
            }
            unit_rotation(cosine, sine)
        })
        .collect()
}

fn unit_rotation(real: f64, imaginary: f64) -> (f64, f64) {
    let norm = real.hypot(imaginary);
    if norm > 1e-12 {
        (real / norm, imaginary / norm)
    } else {
        (1., 0.)
    }
}

fn rotate(point: Point, rotation: (f64, f64)) -> Point {
    Point {
        x: rotation.0 * point.x - rotation.1 * point.y,
        y: rotation.1 * point.x + rotation.0 * point.y,
    }
}

// Rotation smoothness is defined on material edges, not between control pins.
// Area scaling makes the same shape behave identically at another resolution.
fn component_areas(
    original: &[Point],
    triangles: &[[usize; 3]],
    adjacency: &[Vec<(usize, f64)>],
) -> Vec<f64> {
    let mut components = vec![usize::MAX; original.len()];
    let mut areas = Vec::<f64>::new();
    for start in 0..original.len() {
        if components[start] != usize::MAX {
            continue;
        }
        let id = areas.len();
        areas.push(0.);
        let mut queue = vec![start];
        components[start] = id;
        while let Some(i) = queue.pop() {
            for &(j, _) in &adjacency[i] {
                if components[j] == usize::MAX {
                    components[j] = id;
                    queue.push(j);
                }
            }
        }
    }
    for &[a, b, c] in triangles {
        let u = original[b].sub(original[a]);
        let v = original[c].sub(original[a]);
        areas[components[a]] += (u.x * v.y - u.y * v.x).abs() * 0.5;
    }
    (0..original.len()).map(|i| areas[components[i]]).collect()
}

// Penalize jumps of the full material derivative across shared mesh edges.
// This is a quadratic mesh energy, not a rigid disk around each position pin.
// Component area makes the coefficient invariant under uniform image scaling.
fn gradient_regularization(
    original: &[Point],
    triangles: &[[usize; 3]],
    base: &[Vec<(usize, f64)>],
    combined: &[Vec<(usize, f64)>],
    strength: f64,
) -> Vec<Vec<(usize, f64)>> {
    let valid = elements(original, triangles);
    let triangles = valid.iter().map(|e| e.vertices).collect::<Vec<_>>();
    let area = component_areas(original, &triangles, base);
    let gradients = valid
        .iter()
        .map(|e| {
            [
                (e.vertices[0], e.gradients[0]),
                (e.vertices[1], e.gradients[1]),
                (e.vertices[2], e.gradients[2]),
            ]
        })
        .collect::<Vec<_>>();
    let mut edges = combined
        .iter()
        .map(|row| row.iter().copied().collect::<BTreeMap<_, _>>())
        .collect::<Vec<_>>();
    let mut incident = BTreeMap::<(usize, usize), usize>::new();
    for (t, &[a, b, c]) in triangles.iter().enumerate() {
        for (i, j) in [(a, b), (b, c), (c, a)] {
            let key = (i.min(j), i.max(j));
            if let Some(&other) = incident.get(&key) {
                let mut difference = BTreeMap::<usize, Point>::new();
                for &(k, g) in &gradients[t] {
                    difference.entry(k).or_default().x += g.x;
                    difference.entry(k).or_default().y += g.y;
                }
                for &(k, g) in &gradients[other] {
                    difference.entry(k).or_default().x -= g.x;
                    difference.entry(k).or_default().y -= g.y;
                }
                let weight = strength * area[i];
                for (&p, &gp) in &difference {
                    for (&q, &gq) in &difference {
                        if p != q {
                            *edges[p].entry(q).or_default() -= weight * (gp.x * gq.x + gp.y * gq.y);
                        }
                    }
                }
            } else {
                incident.insert(key, t);
            }
        }
    }
    edges
        .into_iter()
        .map(|row| row.into_iter().collect())
        .collect()
}

struct Element {
    vertices: [usize; 3],
    gradients: [Point; 3],
    area: f64,
    transverse: Point,
}

// Minimize area-weighted Cauchy-Riemann residuals to initialize the rotation field.
// All coefficients depend on the rest mesh; changing targets cannot switch
// between nonlinear local minima. Explicit pose handles enter the same solve.
#[allow(clippy::too_many_arguments)]
fn conformal_pose(
    original: &[Point],
    elements: &[Element],
    adjacency: &[Vec<(usize, f64)>],
    handles: &[HandleCell],
    constraints: &HashMap<usize, Point>,
    initial: &[Point],
    widths: &[f64],
) -> Vec<Point> {
    let count = original.len();
    let mut edges = vec![BTreeMap::<usize, f64>::new(); count * 2];
    for e in elements {
        let mut row = Vec::new();
        for k in 0..3 {
            let i = e.vertices[k];
            let g = e.gradients[k];
            row.push((2 * i, g.x, g.y));
            row.push((2 * i + 1, -g.y, g.x));
        }
        for &(i, a, b) in &row {
            for &(j, c, d) in &row {
                if i != j {
                    *edges[i].entry(j).or_default() -= e.area * (a * c + b * d);
                }
            }
        }
    }
    let mut distance = vec![f64::INFINITY; count];
    let mut queue = BinaryHeap::new();
    for &i in constraints.keys() {
        distance[i] = 0.;
        queue.push(Reverse((0_u64, i)));
    }
    while let Some(Reverse((bits, i))) = queue.pop() {
        let d = f64::from_bits(bits);
        if d > distance[i] {
            continue;
        }
        for &(j, _) in &adjacency[i] {
            let next = d + squared_distance(original[i], original[j]).sqrt();
            if next < distance[j] {
                distance[j] = next;
                queue.push(Reverse((next.to_bits(), j)));
            }
        }
    }
    let explicit = handles.iter().any(|h| h.angle.is_some());
    let bias = distance
        .iter()
        .zip(widths)
        .map(|(&d, &width)| {
            if explicit {
                1e-6
            } else {
                1e-6 + 0.05 * (d / width.max(1e-12) - 0.75).clamp(0., 1.).powi(2)
            }
        })
        .collect::<Vec<_>>();
    let mut extra = vec![Point::default(); count * 2];
    for i in 0..count {
        // Resolve the similarity nullspace of a singly pinned component without
        // a rotation constraint. This tiny quadratic term selects its rest pose.
        for &(j, w) in &adjacency[i] {
            let w = w * (bias[i] + bias[j]) * 0.5;
            for axis in 0..2 {
                *edges[2 * i + axis].entry(2 * j + axis).or_default() += w;
            }
            extra[2 * i].x += w * (initial[i].x - initial[j].x);
            extra[2 * i + 1].x += w * (initial[i].y - initial[j].y);
        }
    }
    for h in handles.iter().filter(|h| h.angle.is_some()) {
        let rotation = h.rotation(original, initial);
        for &(j, w) in &h.samples {
            let edge = rotate(original[j].sub(original[h.center]), rotation).scale(h.scale * w);
            for axis in 0..2 {
                *edges[2 * j + axis].entry(2 * h.center + axis).or_default() += w;
                *edges[2 * h.center + axis].entry(2 * j + axis).or_default() += w;
            }
            extra[2 * j].x += edge.x;
            extra[2 * j + 1].x += edge.y;
            extra[2 * h.center].x -= edge.x;
            extra[2 * h.center + 1].x -= edge.y;
        }
    }
    let graph = edges
        .into_iter()
        .map(|row| row.into_iter().collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut fixed = HashMap::new();
    for (&i, &q) in constraints {
        fixed.insert(2 * i, Point { x: q.x, y: 0. });
        fixed.insert(2 * i + 1, Point { x: q.y, y: 0. });
    }
    let first = initial
        .iter()
        .flat_map(|p| [Point { x: p.x, y: 0. }, Point { x: p.y, y: 0. }])
        .collect::<Vec<_>>();
    let result = GlobalSystem::new(&graph, &fixed).step_extra(
        &first,
        &vec![Vec::new(); count * 2],
        &vec![(1., 0.); count * 2],
        &fixed,
        &first,
        &extra,
    );
    result
        .chunks_exact(2)
        .map(|p| Point {
            x: p[0].x,
            y: p[1].x,
        })
        .collect()
}
fn elements(original: &[Point], triangles: &[[usize; 3]]) -> Vec<Element> {
    triangles
        .iter()
        .filter_map(|&[a, b, c]| {
            let u = original[b].sub(original[a]);
            let v = original[c].sub(original[a]);
            let det = u.x * v.y - u.y * v.x;
            if det.abs() < 1e-14 {
                return None;
            }
            Some(Element {
                vertices: [a, b, c],
                area: det.abs() * 0.5,
                transverse: Point { x: 1.0, y: 0.0 },
                gradients: [
                    Point {
                        x: u.y - v.y,
                        y: v.x - u.x,
                    }
                    .scale(1. / det),
                    Point { x: v.y, y: -v.x }.scale(1. / det),
                    Point { x: -u.y, y: u.x }.scale(1. / det),
                ],
            })
        })
        .collect()
}
fn element_adjacency(elements: &[Element], count: usize) -> Vec<Vec<(usize, f64)>> {
    let mut edges = vec![BTreeMap::<usize, f64>::new(); count];
    for e in elements {
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    let a = e.gradients[i];
                    let b = e.gradients[j];
                    *edges[e.vertices[i]].entry(e.vertices[j]).or_default() -=
                        e.area * (a.x * b.x + a.y * b.y);
                }
            }
        }
    }
    edges
        .into_iter()
        .map(|row| row.into_iter().collect())
        .collect()
}
fn material_similarities(
    elements: &[Element],
    pose: &[Point],
    smoothing: &SimilaritySmoother,
) -> Vec<Point> {
    let values = elements
        .iter()
        .map(|e| {
            let mut z = Point::default();
            let (mut fx, mut fy) = (Point::default(), Point::default());
            for i in 0..3 {
                let q = pose[e.vertices[i]].sub(pose[e.vertices[0]]);
                let g = e.gradients[i];
                fx = fx.add(q.scale(g.x));
                fy = fy.add(q.scale(g.y));
                z.x += (q.x * g.x + q.y * g.y) * 0.5;
                z.y += (q.y * g.x - q.x * g.y) * 0.5;
            }
            let determinant = fx.x * fy.y - fx.y * fy.x;
            let norm = fx.x * fx.x + fx.y * fx.y + fy.x * fy.x + fy.y * fy.y;
            z.scale((2. * determinant / norm.max(1e-12)).clamp(0., 1.))
        })
        .collect::<Vec<_>>();
    smoothing.apply(&values)
}
fn element_rhs(
    elements: &[Element],
    deformed: &[Point],
    similarities: &[Point],
    scales: &[f64],
    strain_weight: f64,
    protected: &[bool],
    thickness_weights: &[f64],
) -> Vec<Point> {
    let mut rhs = vec![Point::default(); deformed.len()];
    let mut derivatives = Vec::new();
    for e in elements {
        let (mut a, mut b, mut c, mut d) = (0., 0., 0., 0.);
        for i in 0..3 {
            let q = deformed[e.vertices[i]].sub(deformed[e.vertices[0]]);
            let g = e.gradients[i];
            a += q.x * g.x;
            b += q.x * g.y;
            c += q.y * g.x;
            d += q.y * g.y;
        }
        derivatives.push((a, b, c, d));
    }
    for (index, ((e, z), f)) in elements
        .iter()
        .zip(similarities)
        .zip(derivatives)
        .enumerate()
    {
        let norm = z.x.hypot(z.y);
        let rotation = if norm > 1e-12 {
            (z.x / norm, z.y / norm)
        } else {
            (1., 0.)
        };
        let scale = e.vertices.iter().map(|&i| scales[i]).sum::<f64>() / 3.;
        let (a, b, c, d) = f;
        let (mut tx, mut ty) = (Point::default(), Point::default());
        // A pin must not select a reflected local frame. Its incident elements
        // use the coherent material rotation found before length restoration.
        if protected[index] {
            tx = rotate(Point { x: scale, y: 0. }, rotation);
            ty = rotate(Point { x: 0., y: scale }, rotation);
        } else {
            let angle = 0.5 * (2. * (a * b + c * d)).atan2(a * a + c * c - b * b - d * d);
            let axis = Point {
                x: angle.cos(),
                y: angle.sin(),
            };
            for v in [
                axis,
                Point {
                    x: -axis.y,
                    y: axis.x,
                },
            ] {
                let q = Point {
                    x: a * v.x + b * v.y,
                    y: c * v.x + d * v.y,
                };
                let length = q.x.hypot(q.y);
                // F / max(sigma, epsilon) is continuous through a fold. Normal
                // stretches retain exact unit/explicit scale; a vanishing axis
                // does not jump between opposite unit vectors.
                let target = q.scale(scale / length.max(POLAR_TRANSITION * scale).max(1e-12));
                tx = tx.add(target.scale(v.x));
                ty = ty.add(target.scale(v.y));
            }
            // Restore transverse thickness without inflating a shortened limb
            // to compensate for its lost area. The rest cross section comes
            // from the silhouette, independently of pins and their targets.
            let n = e.transverse;
            let tangent = Point { x: n.y, y: -n.x };
            let longitudinal = Point {
                x: a * tangent.x + b * tangent.y,
                y: c * tangent.x + d * tangent.y,
            };
            let transverse = Point {
                x: a * n.x + b * n.y,
                y: c * n.x + d * n.y,
            };
            let length = longitudinal
                .x
                .hypot(longitudinal.y)
                .max(0.1 * scale)
                .max(1e-12);
            let normal = Point {
                x: -longitudinal.y / length,
                y: longitudinal.x / length,
            };
            let width = normal.x * transverse.x + normal.y * transverse.y;
            let orientation = width / width.abs().max(0.1 * scale).max(1e-12);
            let axial_gradient = Point {
                x: transverse.y / length,
                y: -transverse.x / length,
            }
            .sub(longitudinal.scale(width / (length * length)));
            let gx = axial_gradient.scale(tangent.x).add(normal.scale(n.x));
            let gy = axial_gradient.scale(tangent.y).add(normal.scale(n.y));
            let norm = gx.x * gx.x + gx.y * gx.y + gy.x * gy.x + gy.y * gy.y;
            let pressure = thickness_weights[index] * (scale * orientation - width) / norm.max(1.0);
            tx = tx.add(gx.scale(pressure));
            ty = ty.add(gy.scale(pressure));
        }
        for i in 0..3 {
            let j = e.vertices[i];
            let g = e.gradients[i];
            let target = rotate(g, rotation)
                .scale(scale)
                .add(tx.scale(g.x).add(ty.scale(g.y)).scale(strain_weight));
            rhs[j] = rhs[j].add(target.scale(e.area));
        }
    }
    rhs
}

// Average complex similarities before normalizing them to SO(2). A nearly
// collapsed triangle then carries little directional confidence, instead of
// independently choosing a noisy unit rotation. Use mesh adjacency only.
struct SimilaritySmoother {
    system: GlobalSystem,
    areas: Vec<f64>,
}
impl SimilaritySmoother {
    fn new(elements: &[Element], original: &[Point]) -> Self {
        let n = elements.len();
        let mut graph = vec![BTreeMap::<usize, f64>::new(); n + 1];
        let triangles = elements.iter().map(|e| e.vertices).collect::<Vec<_>>();
        let area = component_areas(
            original,
            &triangles,
            &element_adjacency(elements, original.len()),
        );
        let mut incident = BTreeMap::new();
        let centers = elements
            .iter()
            .map(|e| {
                e.vertices
                    .iter()
                    .fold(Point::default(), |p, &i| p.add(original[i]))
                    .scale(1. / 3.)
            })
            .collect::<Vec<_>>();
        for (t, e) in elements.iter().enumerate() {
            graph[t].insert(n, e.area);
            graph[n].insert(t, e.area);
            let [a, b, c] = e.vertices;
            for (i, j) in [(a, b), (b, c), (c, a)] {
                let key = (i.min(j), i.max(j));
                if let Some(&other) = incident.get(&key) {
                    let w = ROTATION_DIFFUSION
                        * area[i]
                        * squared_distance(original[i], original[j]).sqrt()
                        / squared_distance(centers[t], centers[other])
                            .sqrt()
                            .max(1e-12);
                    graph[t].insert(other, w);
                    graph[other].insert(t, w);
                } else {
                    incident.insert(key, t);
                }
            }
        }
        let graph = graph
            .into_iter()
            .map(|row| row.into_iter().collect())
            .collect::<Vec<_>>();
        let fixed = HashMap::from([(n, Point::default())]);
        Self {
            system: GlobalSystem::new(&graph, &fixed),
            areas: elements.iter().map(|e| e.area).collect(),
        }
    }
    fn apply(&self, values: &[Point]) -> Vec<Point> {
        let n = self.areas.len();
        let mut extra = values
            .iter()
            .zip(&self.areas)
            .map(|(&z, &area)| z.scale(area))
            .collect::<Vec<_>>();
        extra.push(Point::default());
        // The dummy node is fixed at zero, so the homogeneous material response
        // is exactly the screened diffusion solve. Reuse its factor each iteration.
        let mut result = self.system.response(&extra);
        result.truncate(n);
        result
    }
}

// The constrained matrix is constant for all local/global iterations. Build it
// once, with compact indices; small systems reuse a Cholesky factor and larger
// systems use a cached sparse Cholesky factor. A bounded IC-CG fallback is
// reserved for a failed or excessive-fill factorization, independent of targets.
struct GlobalSystem {
    free: Vec<usize>,
    fixed: Vec<Point>,
    solver: PreparedSystem,
}
impl GlobalSystem {
    fn new(adjacency: &[Vec<(usize, f64)>], constraints: &HashMap<usize, Point>) -> Self {
        let free = (0..adjacency.len())
            .filter(|i| !constraints.contains_key(i))
            .collect::<Vec<_>>();
        let mut map = vec![usize::MAX; adjacency.len()];
        for (row, &i) in free.iter().enumerate() {
            map[i] = row;
        }
        let mut rows = Vec::new();
        let mut diagonal = Vec::new();
        let mut fixed = Vec::new();
        for &i in &free {
            let mut row = Vec::new();
            let mut d = 0.0;
            let mut target = Point::default();
            for &(j, w) in &adjacency[i] {
                d += w;
                if map[j] != usize::MAX {
                    row.push((map[j], w));
                } else if let Some(&q) = constraints.get(&j) {
                    target = target.add(q.scale(w));
                }
            }
            rows.push(row);
            diagonal.push(d.max(1.0e-20));
            fixed.push(target);
        }
        let solver = PreparedSystem::new(rows, diagonal, SolveOptions::default());
        Self {
            free,
            fixed,
            solver,
        }
    }
    // Homogeneous boundary response in the same material stiffness metric.
    // Unlike moving just a triangle's vertices, this spreads an area correction
    // through the surface without introducing an isolated fold at its boundary.
    fn response(&self, force: &[Point]) -> Vec<Point> {
        #[cfg(test)]
        WORK.with(|w| w.borrow_mut()[0] += 1);
        let solve = |rhs: &[f64]| self.solver.solve(rhs, vec![0.; self.free.len()]);
        let x = solve(&self.free.iter().map(|&i| force[i].x).collect::<Vec<_>>());
        let y = solve(&self.free.iter().map(|&i| force[i].y).collect::<Vec<_>>());
        let mut result = vec![Point::default(); force.len()];
        for (r, &i) in self.free.iter().enumerate() {
            result[i] = Point { x: x[r], y: y[r] };
        }
        result
    }
    #[cfg(test)]
    fn step(
        &self,
        original: &[Point],
        adjacency: &[Vec<(usize, f64)>],
        rotations: &[(f64, f64)],
        constraints: &HashMap<usize, Point>,
        previous: &[Point],
    ) -> Vec<Point> {
        self.step_extra(
            original,
            adjacency,
            rotations,
            constraints,
            previous,
            &vec![Point::default(); original.len()],
        )
    }
    #[allow(clippy::too_many_arguments)]
    fn step_extra(
        &self,
        original: &[Point],
        adjacency: &[Vec<(usize, f64)>],
        rotations: &[(f64, f64)],
        constraints: &HashMap<usize, Point>,
        previous: &[Point],
        extra: &[Point],
    ) -> Vec<Point> {
        #[cfg(test)]
        WORK.with(|w| w.borrow_mut()[0] += 1);
        let mut bx = vec![0.0; self.free.len()];
        let mut by = bx.clone();
        for (row, &i) in self.free.iter().enumerate() {
            let mut rhs = self.fixed[row].add(extra[i]);
            for &(j, w) in &adjacency[i] {
                let edge = original[i].sub(original[j]);
                rhs = rhs.add(
                    rotate(edge, rotations[i])
                        .add(rotate(edge, rotations[j]))
                        .scale(0.5 * w),
                );
            }
            bx[row] = rhs.x;
            by[row] = rhs.y;
        }
        let solve = |rhs: &[f64], initial: Vec<f64>| self.solver.solve(rhs, initial);
        let x = solve(&bx, self.free.iter().map(|&i| previous[i].x).collect());
        let y = solve(&by, self.free.iter().map(|&i| previous[i].y).collect());
        let mut result = previous.to_vec();
        for (row, &i) in self.free.iter().enumerate() {
            result[i] = Point {
                x: x[row],
                y: y[row],
            };
        }
        for (&i, &q) in constraints {
            result[i] = q;
        }
        result
    }
}
#[cfg(test)]
mod tests {
    #[test]
    fn rotation_diffusion_is_continuous_when_neighboring_frames_are_opposite() {
        use super::*;
        let original = points(&[0., 0., 10., 0., 10., 10., 0., 10.]).unwrap();
        let elements = elements(&original, &[[0, 1, 2], [0, 2, 3]]);
        let smoother = SimilaritySmoother::new(&elements, &original);
        let sample = |offset: f64| {
            let angle = std::f64::consts::PI + offset;
            smoother.apply(&[
                Point { x: 1., y: 0. },
                Point {
                    x: angle.cos(),
                    y: angle.sin(),
                },
            ])
        };
        let before = sample(-1e-8);
        let after = sample(1e-8);
        for (a, b) in before.iter().zip(after) {
            assert!(
                squared_distance(*a, b) < 1e-12,
                "opposite frames introduced a branch jump: {a:?} / {b:?}"
            );
        }
    }

    #[test]
    fn rotation_diffusion_is_equivariant_across_multiple_full_turns() {
        use super::*;
        let original = points(&[0., 0., 10., 0., 10., 10., 0., 10.]).unwrap();
        let elements = elements(&original, &[[0, 1, 2], [0, 2, 3]]);
        let smoother = SimilaritySmoother::new(&elements, &original);
        let values = [
            Point {
                x: 0.2_f64.cos(),
                y: 0.2_f64.sin(),
            },
            Point {
                x: (-0.2_f64).cos(),
                y: (-0.2_f64).sin(),
            },
        ];
        let reference = smoother.apply(&values);
        for degrees in -1080..=1080 {
            let angle = (degrees as f64).to_radians();
            let rotation = (angle.cos(), angle.sin());
            let turned = values.map(|p| rotate(p, rotation));
            let output = smoother.apply(&turned);
            for (a, b) in output.iter().zip(&reference) {
                assert!(
                    squared_distance(*a, rotate(*b, rotation)) < 1e-20,
                    "angle {degrees}"
                );
            }
        }
    }

    #[test]
    fn free_material_rotates_about_its_joint_without_centroid_drift() {
        use super::*;
        let original = vec![
            Point { x: -1., y: -3. },
            Point { x: 1., y: -3. },
            Point { x: 1., y: -1. },
            Point { x: -1., y: -1. },
            Point { x: 0., y: 0. },
            Point { x: -4., y: 1. },
            Point { x: 4., y: 1. },
            Point { x: 0., y: 5. },
        ];
        let part = HandleCell {
            frame: None,
            center: 0,
            samples: vec![(1, 1.), (2, 1.), (3, 1.)],
            angle: None,
            scale: 1.,
        };
        for degrees in [0.0_f64, 83., -170.] {
            let r = (degrees.to_radians().cos(), degrees.to_radians().sin());
            let translation = Point { x: 7., y: 2. };
            for lift in [-2., 0., 2.] {
                let constraints = (4..8)
                    .map(|i| {
                        let mut p = original[i];
                        if i == 5 || i == 6 {
                            p.y += lift;
                        }
                        (i, rotate(p, r).add(translation))
                    })
                    .collect::<HashMap<_, _>>();
                for bias in [-0.3_f64, 0.4] {
                    let mut output = original
                        .iter()
                        .map(|&p| rotate(rotate(p, (bias.cos(), bias.sin())), r).add(translation))
                        .collect::<Vec<_>>();
                    for (&i, &q) in &constraints {
                        output[i] = q;
                    }
                    project_material_parts(
                        &original,
                        &[[0, 1, 2], [0, 2, 3], [3, 2, 4]],
                        &constraints,
                        std::slice::from_ref(&part),
                        &mut output,
                    );
                    let edge = output[1].sub(output[0]);
                    assert!(squared_distance(edge, rotate(Point { x: 2., y: 0. }, r)) < 1e-20);
                    for i in 0..4 {
                        assert!(
                            squared_distance(output[i], rotate(original[i], r).add(translation))
                                < 1e-20,
                            "part translated independently of its joint"
                        );
                    }
                    for (&i, &q) in &constraints {
                        assert_eq!(output[i], q);
                    }
                }
            }
        }
    }

    #[test]
    fn material_projection_is_rotation_equivariant_and_independent_of_diagonals() {
        use super::*;
        let original = vec![
            Point { x: -2., y: -1. },
            Point { x: 2., y: -1. },
            Point { x: 2., y: 1. },
            Point { x: -2., y: 1. },
        ];
        let part = HandleCell {
            frame: None,
            center: 0,
            samples: vec![(1, 1.), (2, 1.), (3, 1.)],
            angle: None,
            scale: 1.,
        };
        for pinned in [false, true] {
            let mut baseline: Option<Vec<Point>> = None;
            for angle in [0.0_f64, -170., -65., 82., 177.] {
                let r = (angle.to_radians().cos(), angle.to_radians().sin());
                let translation = Point { x: 7., y: -3. };
                let field = original
                    .iter()
                    .map(|p| {
                        rotate(
                            Point {
                                x: 1.3 * p.x + 0.4 * p.y,
                                y: 0.8 * p.y,
                            },
                            r,
                        )
                        .add(translation)
                    })
                    .collect::<Vec<_>>();
                let fixed = if pinned {
                    HashMap::from([(0, field[0])])
                } else {
                    HashMap::new()
                };
                let mut results = Vec::new();
                for triangles in [vec![[0, 1, 2], [0, 2, 3]], vec![[0, 1, 3], [1, 2, 3]]] {
                    let mut output = field.clone();
                    project_material_parts(
                        &original,
                        &triangles,
                        &fixed,
                        std::slice::from_ref(&part),
                        &mut output,
                    );
                    if pinned {
                        assert_eq!(output[0], field[0]);
                    }
                    results.push(
                        output
                            .iter()
                            .map(|p| rotate(p.sub(translation), (r.0, -r.1)))
                            .collect::<Vec<_>>(),
                    );
                }
                for (a, b) in results[0].iter().zip(&results[1]) {
                    assert!(
                        squared_distance(*a, *b) < 1e-20,
                        "triangulation biased material rotation"
                    );
                }
                if let Some(reference) = &baseline {
                    for (a, b) in results[0].iter().zip(reference) {
                        assert!(
                            squared_distance(*a, *b) < 1e-20,
                            "world direction biased material rotation"
                        );
                    }
                } else {
                    baseline = Some(results[0].clone());
                }
            }
        }
    }
    use super::*;

    #[test]
    fn finite_handle_local_global_steps_minimize_the_declared_energy() {
        let p = points(&[0., 0., 10., 0., 10., 10., 0., 10., 4., 5.]).unwrap();
        let triangles = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];
        let base = adjacency(&p, &triangles);
        let poses =
            crate::pose::parse(&[0., 0., -8., -3., 0., 1., 0.4, 1.2, 8., 1., 1., 1.]).unwrap();
        let cells = handle_cells(&p, &triangles, &base, &poses);
        let constraints = HashMap::from([(0, Point { x: -8., y: -3. }), (2, p[2])]);
        let system = GlobalSystem::new(&combined_adjacency(&base, &cells), &constraints);
        let energy = |q: &[Point]| {
            let rotations = local_step(&p, q, &base);
            let surface = base
                .iter()
                .enumerate()
                .map(|(i, row)| {
                    row.iter()
                        .map(|&(j, w)| {
                            0.5 * w
                                * squared_distance(
                                    q[i].sub(q[j]),
                                    rotate(p[i].sub(p[j]), rotations[i]),
                                )
                        })
                        .sum::<f64>()
                })
                .sum::<f64>();
            surface
                + cells
                    .iter()
                    .map(|h| {
                        h.samples
                            .iter()
                            .map(|&(j, w)| {
                                w * squared_distance(
                                    q[j].sub(q[h.center]),
                                    rotate(p[j].sub(p[h.center]), h.rotation(&p, q)).scale(h.scale),
                                )
                            })
                            .sum::<f64>()
                    })
                    .sum::<f64>()
        };
        let mut q = p.clone();
        for (&i, &target) in &constraints {
            q[i] = target;
        }
        let initial = energy(&q);
        let mut previous = initial;
        for _ in 0..64 {
            let mut extra = vec![Point::default(); p.len()];
            for h in &cells {
                for &(j, w) in &h.samples {
                    let edge = rotate(p[j].sub(p[h.center]), h.rotation(&p, &q)).scale(h.scale * w);
                    extra[j] = extra[j].add(edge);
                    extra[h.center] = extra[h.center].sub(edge);
                }
            }
            q = system.step_extra(
                &p,
                &base,
                &local_step(&p, &q, &base),
                &constraints,
                &q,
                &extra,
            );
            let next = energy(&q);
            assert!(next <= previous + 1e-8);
            previous = next;
            for (&i, &target) in &constraints {
                assert_eq!(q[i], target);
            }
        }
        assert!(previous < initial * 0.9);
    }

    #[test]
    fn large_sparse_global_system_matches_linear_solution() {
        let n = 420;
        let original = (0..n)
            .map(|i| Point {
                x: i as f64,
                y: 0.0,
            })
            .collect::<Vec<_>>();
        let mut adjacency = vec![Vec::new(); n];
        for i in 1..n {
            adjacency[i - 1].push((i, 1.0));
            adjacency[i].push((i - 1, 1.0));
        }
        let constraints = HashMap::from([
            (0, Point { x: 0.0, y: 0.0 }),
            (
                n - 1,
                Point {
                    x: (n - 1) as f64 * 2.0,
                    y: 10.0,
                },
            ),
        ]);
        let system = GlobalSystem::new(&adjacency, &constraints);
        assert!(system.solver.uses_sparse_factor());
        let result = system.step(
            &original,
            &adjacency,
            &vec![(1.0, 0.0); n],
            &constraints,
            &original,
        );
        for (i, p) in result.iter().enumerate() {
            assert!((p.x - i as f64 * 2.0).abs() < 1.0e-5);
            assert!((p.y - 10.0 * i as f64 / (n - 1) as f64).abs() < 1.0e-5);
        }
    }

    #[test]
    fn local_global_steps_decrease_rigidity_energy() {
        let original = points(&[0., 0., 10., 0., 10., 10., 0., 10., 4., 5.]).unwrap();
        let adjacency = adjacency(&original, &[[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]);
        let constraints = HashMap::from([
            (0, original[0]),
            (1, Point { x: 13., y: 3. }),
            (2, Point { x: 8., y: 15. }),
        ]);
        let system = GlobalSystem::new(&adjacency, &constraints);
        let mut q = original.clone();
        for (&i, &p) in &constraints {
            q[i] = p;
        }
        let energy = |q: &[Point]| {
            let rotations = local_step(&original, q, &adjacency);
            adjacency
                .iter()
                .enumerate()
                .map(|(i, edges)| {
                    edges
                        .iter()
                        .map(|&(j, w)| {
                            w * squared_distance(
                                q[i].sub(q[j]),
                                rotate(original[i].sub(original[j]), rotations[i]),
                            )
                        })
                        .sum::<f64>()
                })
                .sum::<f64>()
        };
        let initial = energy(&q);
        let mut previous = initial;
        for _ in 0..24 {
            q = system.step(
                &original,
                &adjacency,
                &local_step(&original, &q, &adjacency),
                &constraints,
                &q,
            );
            let next = energy(&q);
            assert!(next <= previous + 1.0e-8);
            previous = next;
            for (&i, &p) in &constraints {
                assert_eq!(q[i], p);
            }
        }
        assert!(previous < initial * 0.9);
    }

    #[test]
    fn hard_constraints_are_exact() {
        let vertices = [0.0, 0.0, 10.0, 0.0, 10.0, 10.0, 0.0, 10.0];
        let indices = [0, 1, 2, 0, 2, 3];
        let result = solve(
            &vertices,
            &indices,
            &[0.0, 0.0, 10.0, 10.0],
            &[2.0, 3.0, 13.0, 14.0],
        )
        .unwrap();
        assert_eq!(result[0], Point { x: 2.0, y: 3.0 });
        assert_eq!(result[2], Point { x: 13.0, y: 14.0 });
        assert!(
            result
                .iter()
                .all(|point| point.x.is_finite() && point.y.is_finite())
        );
    }

    #[test]
    fn unchanged_constraints_preserve_a_triangle() {
        let vertices = [0.0, 0.0, 10.0, 0.0, 0.0, 10.0];
        let result = solve(&vertices, &[0, 1, 2], &vertices, &vertices).unwrap();
        for (actual, expected) in result.iter().zip(points(&vertices).unwrap()) {
            assert!(squared_distance(*actual, expected) < 1.0e-12);
        }
    }
}
