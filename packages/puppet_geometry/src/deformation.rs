use std::{cmp::Reverse, collections::BinaryHeap};

use aviutl2::anyhow::{self, Context as _};

const EPSILON: f64 = 1.0e-10;

use crate::{math::Point, mls::RigidMls};

pub(crate) enum Method {
    RigidMls(f64),
    Arap,
}

#[derive(Clone, Copy)]
struct Sample {
    original: Point,
    deformed: Point,
    layer: f64,
}

pub(crate) struct Deformation {
    pub vertices: Vec<f64>,
    /// Triangle-list vertices, four values per vertex: x, y, u, v.
    pub render_vertices: Vec<f64>,
    /// Original mesh edges sampled through the same field as the renderer.
    pub wire_vertices: Vec<f64>,
}

#[allow(clippy::too_many_arguments)] // Mirrors the flat AviUtl2 module ABI.
#[cfg(test)]
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
    deform(
        vertices,
        indices,
        pin_source,
        pin_destination,
        pin_layers,
        divisions,
        width,
        height,
        &[],
        Method::RigidMls(stiffness),
    )
}

#[allow(clippy::too_many_arguments)] // Mirrors the flat AviUtl2 module ABI.
#[cfg(test)]
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
    deform(
        vertices,
        indices,
        pin_source,
        pin_destination,
        pin_layers,
        divisions,
        width,
        height,
        &[],
        Method::Arap,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn deform(
    vertices: &[f64],
    indices: &[i32],
    pin_source: &[f64],
    pin_destination: &[f64],
    pin_layers: &[f64],
    divisions: i32,
    width: f64,
    height: f64,
    pose_data: &[f64],
    method: Method,
) -> anyhow::Result<Deformation> {
    if let Method::RigidMls(exponent) = method {
        anyhow::ensure!(exponent.is_finite() && exponent >= 0., "invalid stiffness");
    }
    let poses = crate::pose::parse(pose_data)?;
    validate(
        vertices,
        indices,
        pin_source,
        pin_destination,
        pin_layers,
        width,
        height,
    )?;
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
    let original = to_points(vertices);
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
    let overlap_distances = geodesic_distances(&original, &triangles, &overlap_sources);

    let mls = match method {
        Method::RigidMls(exponent) => Some(RigidMls::new(
            &geometry_sources,
            &geometry_destinations,
            &poses,
            exponent,
        )),
        Method::Arap => None,
    };
    let deformed = if let Some(field) = &mls {
        original
            .iter()
            .map(|&p| field.evaluate(p))
            .collect::<Vec<_>>()
    } else {
        crate::arap::solve_with_poses(
            vertices,
            indices,
            &geometry_sources,
            &geometry_destinations,
            &poses,
        )?
    };
    let mut flat_deformed = Vec::with_capacity(vertices.len());
    for point in &deformed {
        flat_deformed.extend([point.x, point.y]);
    }

    if divisions == 0 {
        return Ok(Deformation {
            vertices: flat_deformed,
            render_vertices: Vec::new(),
            wire_vertices: Vec::new(),
        });
    }
    let mut wire_vertices = Vec::new();
    let mut wire_edges = std::collections::HashSet::new();
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
                let deformed = mls.as_ref().map_or_else(
                    || interpolate_points(target, barycentric),
                    |field| field.evaluate(original),
                );
                let sample_distances =
                    interpolate_distances(&overlap_distances, triangle, barycentric);
                sample_row.push(Sample {
                    original,
                    deformed,
                    layer: (original.y + height * 0.5) / height
                        + overlap_layer(&sample_distances, &overlap_layers, &overlap_ranges),
                });
            }
        }
        append_wire(
            &mut wire_vertices,
            &mut wire_edges,
            triangle,
            &samples,
            divisions,
        );
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
        wire_vertices,
        render_vertices: rendered.into_iter().flat_map(|(v, _)| v).collect(),
    })
}

fn to_points(values: &[f64]) -> Vec<Point> {
    values
        .chunks_exact(2)
        .map(|v| Point { x: v[0], y: v[1] })
        .collect()
}

fn validate(
    vertices: &[f64],
    indices: &[i32],
    sources: &[f64],
    destinations: &[f64],
    layers: &[f64],
    width: f64,
    height: f64,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        vertices.len() >= 6 && vertices.len().is_multiple_of(2),
        "invalid vertices"
    );
    anyhow::ensure!(
        indices.len() >= 3 && indices.len().is_multiple_of(3),
        "invalid indices"
    );
    anyhow::ensure!(
        sources.len().is_multiple_of(2) && sources.len() == destinations.len(),
        "invalid pin arrays"
    );
    anyhow::ensure!(
        vertices
            .iter()
            .chain(sources)
            .chain(destinations)
            .chain(layers)
            .all(|v| v.is_finite()),
        "non-finite coordinates"
    );
    anyhow::ensure!(
        width.is_finite() && height.is_finite() && width > 0.0 && height > 0.0,
        "invalid dimensions"
    );
    Ok(())
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

fn interpolate_distances(fields: &[Vec<f64>], triangle: [usize; 3], b: [f64; 3]) -> Vec<f64> {
    fields
        .iter()
        .map(|field| {
            (0..3)
                .filter(|&i| b[i] > 0.0)
                .map(|i| b[i] * field[triangle[i]])
                .sum()
        })
        .collect()
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
            // Nonnegative IEEE-754 distances sort like their unsigned bits.
            let mut queue = BinaryHeap::from([Reverse((initial.to_bits(), start))]);
            while let Some(Reverse((bits, vertex))) = queue.pop() {
                let distance = f64::from_bits(bits);
                if distance > result[vertex] {
                    continue;
                }
                for &(next, edge) in &adjacency[vertex] {
                    let candidate = distance + edge;
                    if candidate < result[next] {
                        result[next] = candidate;
                        queue.push(Reverse((candidate.to_bits(), next)));
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

fn append_wire(
    out: &mut Vec<f64>,
    seen: &mut std::collections::HashSet<(usize, usize)>,
    triangle: [usize; 3],
    samples: &[Vec<Sample>],
    divisions: usize,
) {
    for (edge, (a, b)) in [
        (triangle[0], triangle[1]),
        (triangle[0], triangle[2]),
        (triangle[1], triangle[2]),
    ]
    .into_iter()
    .enumerate()
    {
        if !seen.insert((a.min(b), a.max(b))) {
            continue;
        }
        let sample = |k: usize| match edge {
            0 => samples[0][k].deformed,
            1 => samples[k][0].deformed,
            _ => samples[k][divisions - k].deformed,
        };
        for k in 0..divisions {
            let a = sample(k);
            let b = sample(k + 1);
            out.extend([a.x, a.y, b.x, b.y]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mls_render_samples_evaluate_the_field_without_arap_or_mesh_interpolation() {
        let source = [0., 0., 10., 0., 0., 10.];
        let target = [1., 0., 12., 3., -2., 7.];
        crate::arap::take_solve_calls();
        let output = deform_mls(
            &source,
            &[0, 1, 2],
            &source,
            &target,
            &[0.; 3],
            1.,
            2,
            20.,
            20.,
        )
        .unwrap();
        assert_eq!(crate::arap::take_solve_calls(), 0);
        let q = RigidMls::new(&source, &target, &[], 1.).evaluate(Point { x: 5., y: 0. });
        let sample = output
            .render_vertices
            .chunks_exact(4)
            .find(|p| p[2] == 0.75 && p[3] == 0.5)
            .unwrap();
        assert!((sample[0] - q.x).hypot(sample[1] - q.y) < 1e-12);
        assert!(
            (q.x - 6.5).hypot(q.y - 1.5) > 1e-4,
            "rendering interpolated the vertex solution"
        );
    }

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
    fn overlap_pin_changes_triangle_order_in_both_directions() {
        let vertices = [-10.0, -10.0, 10.0, -10.0, 10.0, 10.0, -10.0, 10.0];
        let indices = [0, 1, 2, 0, 2, 3];
        let bottom_right_is_first = |render: &[f64]| {
            render[..12]
                .chunks_exact(4)
                .any(|v| (v[2] - 1.0).abs() < 1e-8 && v[3].abs() < 1e-8)
        };
        for method in 1..=2 {
            let deform = |layer| {
                if method == 1 {
                    deform_mls(
                        &vertices,
                        &indices,
                        &[10.0, -10.0],
                        &[10.0, -10.0],
                        &[layer, -16.0],
                        1.0,
                        1,
                        20.0,
                        20.0,
                    )
                } else {
                    deform_arap(
                        &vertices,
                        &indices,
                        &[10.0, -10.0],
                        &[10.0, -10.0],
                        &[layer, -16.0],
                        1,
                        20.0,
                        20.0,
                    )
                }
                .unwrap()
            };
            assert!(bottom_right_is_first(&deform(-2.0).render_vertices));
            assert!(!bottom_right_is_first(&deform(2.0).render_vertices));
        }
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
