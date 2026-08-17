use std::{cmp::Ordering, collections::BinaryHeap};

use aviutl2::anyhow::{self, Context as _};

const EPSILON: f64 = 1.0e-10;

#[derive(Clone, Copy, Debug)]
struct Point {
    x: f64,
    y: f64,
}

#[derive(Clone, Copy)]
struct Sample {
    original: Point,
    deformed: Point,
    layer: f64,
}

#[derive(Clone, Copy)]
struct QueueEntry {
    distance: f64,
    vertex: usize,
}

impl Eq for QueueEntry {}
impl PartialEq for QueueEntry {
    fn eq(&self, other: &Self) -> bool {
        self.distance == other.distance && self.vertex == other.vertex
    }
}
impl Ord for QueueEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .distance
            .total_cmp(&self.distance)
            .then_with(|| other.vertex.cmp(&self.vertex))
    }
}
impl PartialOrd for QueueEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub(crate) struct Deformation {
    pub vertices: Vec<f64>,
    /// Triangle-list vertices, four values per vertex: x, y, u, v.
    pub render_vertices: Vec<f64>,
}

#[allow(clippy::too_many_arguments)] // Mirrors the flat AviUtl2 module ABI.
pub(crate) fn deform_mls(
    vertices: &[f64],
    indices: &[i32],
    pin_source: &[f64],
    pin_destination: &[f64],
    pin_layers: &[f64],
    stiffness: f64,
    divisions: i32,
    width: f64,
    height: f64,
) -> anyhow::Result<Deformation> {
    anyhow::ensure!(
        vertices.len() >= 6 && vertices.len().is_multiple_of(2),
        "invalid vertices"
    );
    anyhow::ensure!(
        indices.len() >= 3 && indices.len().is_multiple_of(3),
        "invalid indices"
    );
    anyhow::ensure!(pin_source.len().is_multiple_of(2), "invalid source pins");
    anyhow::ensure!(
        pin_destination.len() == pin_source.len(),
        "pin arrays differ in length"
    );
    let pin_count = pin_source.len() / 2;
    anyhow::ensure!(
        pin_layers.len() == pin_count || pin_layers.len() == pin_count * 2,
        "invalid pin layers"
    );
    let (pin_layers, pin_ranges) = if pin_layers.len() == pin_count * 2 {
        pin_layers.split_at(pin_count)
    } else {
        (pin_layers, &[][..])
    };
    anyhow::ensure!(
        width > 0.0 && height > 0.0,
        "image dimensions must be positive"
    );
    anyhow::ensure!(
        stiffness.is_finite() && stiffness >= 0.0,
        "invalid stiffness"
    );

    let points = vertices
        .chunks_exact(2)
        .map(|v| Point { x: v[0], y: v[1] })
        .collect::<Vec<_>>();
    let triangles = indices
        .chunks_exact(3)
        .map(|t| {
            let triangle = [
                usize::try_from(t[0]).context("negative mesh index")?,
                usize::try_from(t[1]).context("negative mesh index")?,
                usize::try_from(t[2]).context("negative mesh index")?,
            ];
            anyhow::ensure!(
                triangle.iter().all(|&i| i < points.len()),
                "mesh index out of bounds"
            );
            Ok(triangle)
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let all_sources = to_points(pin_source);
    let all_destinations = to_points(pin_destination);
    let mut sources = Vec::new();
    let mut destinations = Vec::new();
    let mut overlap_sources = Vec::new();
    let mut overlap_layers = Vec::new();
    let mut overlap_ranges = Vec::new();
    for index in 0..pin_count {
        let encoded_range = pin_ranges.get(index).copied().unwrap_or(0.0);
        if encoded_range < 0.0 {
            overlap_sources.push(all_sources[index]);
            overlap_layers.push(pin_layers[index]);
            overlap_ranges.push((-encoded_range - 1.0).max(0.0));
        } else {
            sources.push(all_sources[index]);
            destinations.push(all_destinations[index]);
        }
    }
    if sources.is_empty() {
        sources.push(points[0]);
        destinations.push(points[0]);
    }
    let distances = geodesic_distances(&points, &triangles, &sources);
    let overlap_distances = geodesic_distances(&points, &triangles, &overlap_sources);
    let moved = sources
        .iter()
        .zip(&destinations)
        .any(|(a, b)| a.x != b.x || a.y != b.y);

    let evaluate = |point: Point, distance: &[f64], layer: f64| {
        mls_rigid(
            point,
            distance,
            &sources,
            &destinations,
            stiffness,
            moved,
            layer,
        )
    };
    let deformed = points
        .iter()
        .enumerate()
        .map(|(i, &point)| {
            let d = distances.iter().map(|pin| pin[i]).collect::<Vec<_>>();
            let overlap_d = overlap_distances
                .iter()
                .map(|pin| pin[i])
                .collect::<Vec<_>>();
            evaluate(
                point,
                &d,
                overlap_layer(&overlap_d, &overlap_layers, &overlap_ranges),
            )
        })
        .collect::<Vec<_>>();
    let mut flat_deformed = Vec::with_capacity(points.len() * 2);
    for (point, _) in &deformed {
        flat_deformed.extend([point.x, point.y]);
    }

    let divisions = divisions.clamp(1, 32) as usize;
    let mut rendered = Vec::<([f64; 12], f64)>::new();
    for triangle in triangles {
        let original = [
            points[triangle[0]],
            points[triangle[1]],
            points[triangle[2]],
        ];
        let mut samples = vec![Vec::<Sample>::new(); divisions + 1];
        for (row, sample_row) in samples.iter_mut().enumerate() {
            for col in 0..=divisions - row {
                let b = [
                    1.0 - (row + col) as f64 / divisions as f64,
                    col as f64 / divisions as f64,
                    row as f64 / divisions as f64,
                ];
                let point = interpolate_points(original, b);
                let d = distances
                    .iter()
                    .map(|pin| {
                        b[0] * pin[triangle[0]] + b[1] * pin[triangle[1]] + b[2] * pin[triangle[2]]
                    })
                    .collect::<Vec<_>>();
                let overlap_d = overlap_distances
                    .iter()
                    .map(|pin| {
                        b[0] * pin[triangle[0]] + b[1] * pin[triangle[1]] + b[2] * pin[triangle[2]]
                    })
                    .collect::<Vec<_>>();
                let layer_offset = overlap_layer(&overlap_d, &overlap_layers, &overlap_ranges);
                let (deformed, layer_offset) = evaluate(point, &d, layer_offset);
                sample_row.push(Sample {
                    original: point,
                    deformed,
                    layer: (point.y + height * 0.5) / height + layer_offset,
                });
            }
        }
        for row in 0..divisions {
            for col in 0..divisions - row {
                push_triangle(
                    &mut rendered,
                    [
                        samples[row][col],
                        samples[row][col + 1],
                        samples[row + 1][col],
                    ],
                    width,
                    height,
                );
                if col + 1 < divisions - row {
                    push_triangle(
                        &mut rendered,
                        [
                            samples[row][col + 1],
                            samples[row + 1][col + 1],
                            samples[row + 1][col],
                        ],
                        width,
                        height,
                    );
                }
            }
        }
    }
    rendered.sort_by(|a, b| a.1.total_cmp(&b.1));
    Ok(Deformation {
        vertices: flat_deformed,
        render_vertices: rendered.into_iter().flat_map(|(v, _)| v).collect(),
    })
}

#[allow(clippy::too_many_arguments)] // Mirrors the flat AviUtl2 module ABI.
pub(crate) fn deform_arap(
    vertices: &[f64],
    indices: &[i32],
    pin_source: &[f64],
    pin_destination: &[f64],
    pin_layers: &[f64],
    divisions: i32,
    width: f64,
    height: f64,
) -> anyhow::Result<Deformation> {
    let pin_count = pin_source.len() / 2;
    anyhow::ensure!(
        pin_layers.len() == pin_count || pin_layers.len() == pin_count * 2,
        "invalid pin layers"
    );
    let (pin_layers, pin_ranges) = if pin_layers.len() == pin_count * 2 {
        pin_layers.split_at(pin_count)
    } else {
        (pin_layers, &[][..])
    };
    anyhow::ensure!(
        width > 0.0 && height > 0.0,
        "image dimensions must be positive"
    );
    let original = to_points(vertices);
    anyhow::ensure!(
        original.len() >= 3 && original.len() * 2 == vertices.len(),
        "invalid vertices"
    );
    let triangles = indices
        .chunks_exact(3)
        .map(|t| {
            let triangle = [
                usize::try_from(t[0]).context("negative mesh index")?,
                usize::try_from(t[1]).context("negative mesh index")?,
                usize::try_from(t[2]).context("negative mesh index")?,
            ];
            anyhow::ensure!(
                triangle.iter().all(|&i| i < original.len()),
                "mesh index out of bounds"
            );
            Ok(triangle)
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    anyhow::ensure!(triangles.len() * 3 == indices.len(), "invalid indices");
    let all_sources = to_points(pin_source);
    let mut geometry_sources = Vec::new();
    let mut geometry_destinations = Vec::new();
    let mut overlap_sources = Vec::new();
    let mut overlap_layers = Vec::new();
    let mut overlap_ranges = Vec::new();
    for index in 0..pin_count {
        let encoded_range = pin_ranges.get(index).copied().unwrap_or(0.0);
        if encoded_range < 0.0 {
            overlap_sources.push(all_sources[index]);
            overlap_layers.push(pin_layers[index]);
            overlap_ranges.push((-encoded_range - 1.0).max(0.0));
        } else {
            geometry_sources.extend([pin_source[index * 2], pin_source[index * 2 + 1]]);
            geometry_destinations
                .extend([pin_destination[index * 2], pin_destination[index * 2 + 1]]);
        }
    }
    if geometry_sources.is_empty() {
        geometry_sources.extend([original[0].x, original[0].y]);
        geometry_destinations.extend([original[0].x, original[0].y]);
    }
    let overlap_distances = geodesic_distances(&original, &triangles, &overlap_sources);

    let geometry_points = to_points(&geometry_sources);
    let target_points = to_points(&geometry_destinations);
    let distances = geodesic_distances(&original, &triangles, &geometry_points);
    let moved = geometry_points.iter().zip(&target_points).any(|(a, b)| a.x != b.x || a.y != b.y);

    let mut initial_guess = Vec::with_capacity(original.len());
    for (i, &point) in original.iter().enumerate() {
        let d = distances.iter().map(|pin| pin[i]).collect::<Vec<_>>();
        let (deformed, _) = mls_rigid(point, &d, &geometry_points, &target_points, 1.0, moved, 0.0);
        initial_guess.push(crate::arap::Point { x: deformed.x, y: deformed.y });
    }

    let solved = crate::arap::solve(
        vertices,
        indices,
        &geometry_sources,
        &geometry_destinations,
        &initial_guess,
    )?;
    let deformed = solved
        .iter()
        .map(|p| Point { x: p.x, y: p.y })
        .collect::<Vec<_>>();
    let mut flat_deformed = Vec::with_capacity(vertices.len());
    for point in &deformed {
        flat_deformed.extend([point.x, point.y]);
    }

    let divisions = divisions.clamp(1, 32) as usize;
    let mut rendered = Vec::<([f64; 12], f64)>::new();
    for triangle in triangles {
        let source = [
            original[triangle[0]],
            original[triangle[1]],
            original[triangle[2]],
        ];
        let target = [
            deformed[triangle[0]],
            deformed[triangle[1]],
            deformed[triangle[2]],
        ];
        let mut samples = vec![Vec::<Sample>::new(); divisions + 1];
        for (row, sample_row) in samples.iter_mut().enumerate() {
            for col in 0..=divisions - row {
                let barycentric = [
                    1.0 - (row + col) as f64 / divisions as f64,
                    col as f64 / divisions as f64,
                    row as f64 / divisions as f64,
                ];
                let original = interpolate_points(source, barycentric);
                let deformed = interpolate_points(target, barycentric);
                let sample_distances = overlap_distances
                    .iter()
                    .map(|pin| {
                        barycentric[0] * pin[triangle[0]]
                            + barycentric[1] * pin[triangle[1]]
                            + barycentric[2] * pin[triangle[2]]
                    })
                    .collect::<Vec<_>>();
                sample_row.push(Sample {
                    original,
                    deformed,
                    layer: (original.y + height * 0.5) / height
                        + overlap_layer(&sample_distances, &overlap_layers, &overlap_ranges),
                });
            }
        }
        for row in 0..divisions {
            for col in 0..divisions - row {
                push_triangle(
                    &mut rendered,
                    [
                        samples[row][col],
                        samples[row][col + 1],
                        samples[row + 1][col],
                    ],
                    width,
                    height,
                );
                if col + 1 < divisions - row {
                    push_triangle(
                        &mut rendered,
                        [
                            samples[row][col + 1],
                            samples[row + 1][col + 1],
                            samples[row + 1][col],
                        ],
                        width,
                        height,
                    );
                }
            }
        }
    }
    rendered.sort_by(|a, b| a.1.total_cmp(&b.1));
    Ok(Deformation {
        vertices: flat_deformed,
        render_vertices: rendered.into_iter().flat_map(|(v, _)| v).collect(),
    })
}

fn to_points(values: &[f64]) -> Vec<Point> {
    values
        .chunks_exact(2)
        .map(|v| Point { x: v[0], y: v[1] })
        .collect()
}

fn push_triangle(out: &mut Vec<([f64; 12], f64)>, samples: [Sample; 3], width: f64, height: f64) {
    let mut flat = [0.0; 12];
    for (i, sample) in samples.iter().enumerate() {
        flat[i * 4] = sample.deformed.x;
        flat[i * 4 + 1] = sample.deformed.y;
        flat[i * 4 + 2] = (sample.original.x + width * 0.5) / width;
        flat[i * 4 + 3] = (sample.original.y + height * 0.5) / height;
    }
    let layer = samples.iter().map(|sample| sample.layer).sum::<f64>() / 3.0;
    out.push((flat, layer));
}

fn interpolate_points(points: [Point; 3], b: [f64; 3]) -> Point {
    Point {
        x: b[0] * points[0].x + b[1] * points[1].x + b[2] * points[2].x,
        y: b[0] * points[0].y + b[1] * points[1].y + b[2] * points[2].y,
    }
}

fn geodesic_distances(points: &[Point], triangles: &[[usize; 3]], pins: &[Point]) -> Vec<Vec<f64>> {
    let mut adjacency = vec![Vec::<(usize, f64)>::new(); points.len()];
    for triangle in triangles {
        for (a, b) in [
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ] {
            let distance = distance(points[a], points[b]);
            adjacency[a].push((b, distance));
            adjacency[b].push((a, distance));
        }
    }
    pins.iter()
        .map(|&pin| {
            let (start, initial) = points
                .iter()
                .enumerate()
                .map(|(i, &point)| (i, distance(pin, point)))
                .min_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap();
            let mut result = vec![f64::INFINITY; points.len()];
            result[start] = initial;
            let mut queue = BinaryHeap::from([QueueEntry {
                distance: initial,
                vertex: start,
            }]);
            while let Some(entry) = queue.pop() {
                if entry.distance > result[entry.vertex] {
                    continue;
                }
                for &(next, edge) in &adjacency[entry.vertex] {
                    let candidate = entry.distance + edge;
                    if candidate < result[next] {
                        result[next] = candidate;
                        queue.push(QueueEntry {
                            distance: candidate,
                            vertex: next,
                        });
                    }
                }
            }
            result
        })
        .collect()
}

fn distance(a: Point, b: Point) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

fn mls_rigid(
    point: Point,
    distances: &[f64],
    source: &[Point],
    destination: &[Point],
    stiffness: f64,
    moved: bool,
    layer: f64,
) -> (Point, f64) {
    if !moved {
        return (point, layer);
    }
    if let Some(index) = distances.iter().position(|&d| d < 1.0e-5) {
        return (destination[index], layer);
    }
    let weights = distances
        .iter()
        .map(|d| 1.0 / d.powf(2.0 * stiffness))
        .collect::<Vec<_>>();
    let sum = weights.iter().sum::<f64>();
    if sum < EPSILON || !sum.is_finite() {
        return (point, layer);
    }
    let weighted = |points: &[Point]| Point {
        x: weights
            .iter()
            .zip(points)
            .map(|(w, p)| w * p.x)
            .sum::<f64>()
            / sum,
        y: weights
            .iter()
            .zip(points)
            .map(|(w, p)| w * p.y)
            .sum::<f64>()
            / sum,
    };
    let p_center = weighted(source);
    let q_center = weighted(destination);
    let mut m11 = 0.0;
    let mut m12 = 0.0;
    for ((&weight, p), q) in weights.iter().zip(source).zip(destination) {
        let (px, py) = (p.x - p_center.x, p.y - p_center.y);
        let (qx, qy) = (q.x - q_center.x, q.y - q_center.y);
        m11 += weight * (qx * px + qy * py);
        m12 += weight * (-qx * py + qy * px);
    }
    let norm = m11.hypot(m12);
    if norm < EPSILON {
        return (
            Point {
                x: point.x + q_center.x - p_center.x,
                y: point.y + q_center.y - p_center.y,
            },
            layer,
        );
    }
    let (x, y) = (point.x - p_center.x, point.y - p_center.y);
    (
        Point {
            x: (x * m11 - y * m12) / norm + q_center.x,
            y: (x * m12 + y * m11) / norm + q_center.y,
        },
        layer,
    )
}

fn overlap_layer(distances: &[f64], layers: &[f64], ranges: &[f64]) -> f64 {
    if ranges.is_empty() {
        return 0.0;
    }
    distances
        .iter()
        .zip(layers)
        .zip(ranges)
        .map(|((&distance, &layer), &range)| {
            if range <= EPSILON || layer == 0.0 || distance >= range {
                return 0.0;
            }
            let t = (1.0 - distance / range).clamp(0.0, 1.0);
            layer * t * t * (3.0 - 2.0 * t)
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    const TRIANGLE: [f64; 6] = [0.0, 0.0, 10.0, 0.0, 0.0, 10.0];
    const INDICES: [i32; 3] = [0, 1, 2];

    #[test]
    fn unchanged_pins_preserve_vertices_and_generate_uvs() {
        let pins = [0.0, 0.0, 10.0, 0.0];
        let result = deform_mls(
            &TRIANGLE,
            &INDICES,
            &pins,
            &pins,
            &[0.0, 0.0],
            1.0,
            1,
            20.0,
            20.0,
        )
        .unwrap();
        assert_eq!(result.vertices, TRIANGLE);
        assert_eq!(result.render_vertices.len(), 12);
        assert_eq!(&result.render_vertices[0..4], &[0.0, 0.0, 0.5, 0.5]);
    }

    #[test]
    fn zero_pins_preserve_the_mesh() {
        let mls = deform_mls(&TRIANGLE, &INDICES, &[], &[], &[], 1.0, 1, 20.0, 20.0).unwrap();
        assert_eq!(mls.vertices, TRIANGLE);

        let arap = deform_arap(&TRIANGLE, &INDICES, &[], &[], &[], 1, 20.0, 20.0).unwrap();
        assert_eq!(arap.vertices, TRIANGLE);
    }

    #[test]
    fn subdivision_produces_d_squared_triangles() {
        let pins = [0.0, 0.0];
        let result = deform_mls(
            &TRIANGLE,
            &INDICES,
            &pins,
            &pins,
            &[0.0],
            1.0,
            3,
            20.0,
            20.0,
        )
        .unwrap();
        assert_eq!(result.render_vertices.len(), 3 * 3 * 3 * 4);
    }

    #[test]
    fn pin_snap_moves_its_mesh_vertex_exactly() {
        let result = deform_mls(
            &TRIANGLE,
            &INDICES,
            &[0.0, 0.0],
            &[3.0, 4.0],
            &[0.0],
            1.0,
            1,
            20.0,
            20.0,
        )
        .unwrap();
        assert_eq!(&result.vertices[0..2], &[3.0, 4.0]);
    }

    #[test]
    fn one_mls_pin_translates_the_whole_mesh() {
        let result = deform_mls(
            &TRIANGLE,
            &INDICES,
            &[0.0, 0.0],
            &[3.0, 4.0],
            &[0.0],
            1.0,
            1,
            20.0,
            20.0,
        )
        .unwrap();
        assert_eq!(result.vertices, [3.0, 4.0, 13.0, 4.0, 3.0, 14.0]);
    }

    #[test]
    fn overlap_layer_fades_smoothly_to_zero_at_range() {
        assert_eq!(overlap_layer(&[0.0], &[2.0], &[10.0]), 2.0);
        assert_eq!(overlap_layer(&[5.0], &[2.0], &[10.0]), 1.0);
        assert_eq!(overlap_layer(&[10.0], &[2.0], &[10.0]), 0.0);
    }

    #[test]
    fn overlap_pin_is_not_a_geometry_constraint() {
        let result = deform_mls(
            &TRIANGLE,
            &INDICES,
            &[0.0, 0.0, 10.0, 0.0],
            &[0.0, 0.0, 1_000.0, 1_000.0],
            // First half is layer, second half is encoded range. A negative
            // range marks an overlap-only pin; -11 means an actual range of 10.
            &[0.0, 1.0, 0.0, -11.0],
            1.0,
            1,
            20.0,
            20.0,
        )
        .unwrap();
        assert_eq!(result.vertices, TRIANGLE);
    }

    #[test]
    fn arap_pipeline_returns_complete_subdivided_triangles() {
        let result = deform_arap(
            &TRIANGLE,
            &INDICES,
            &[0.0, 0.0],
            &[2.0, 3.0],
            &[0.0],
            2,
            20.0,
            20.0,
        )
        .unwrap();
        assert_eq!(&result.vertices[0..2], &[2.0, 3.0]);
        assert_eq!(result.render_vertices.len(), 2 * 2 * 3 * 4);
    }
}
