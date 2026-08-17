use std::collections::{HashMap, HashSet};

use aviutl2::anyhow::{self, Context as _};

const ITERATIONS: usize = 8;
const CG_ITERATIONS: usize = 128;
const TOLERANCE: f64 = 1.0e-8;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }

    fn scale(self, value: f64) -> Self {
        Self {
            x: self.x * value,
            y: self.y * value,
        }
    }
}

/// Local/global ARAP with hard positional constraints.
///
/// Uniform edge weights are intentional here: generated puppet meshes can
/// contain very thin boundary triangles, for which cotangent weights become
/// unstable. The solver interface is independent of that choice so weights
/// can be upgraded without changing the Lua/module boundary.
pub(crate) fn solve(
    vertices: &[f64],
    indices: &[i32],
    pin_source: &[f64],
    pin_destination: &[f64],
    initial_guess: &[Point],
) -> anyhow::Result<Vec<Point>> {
    let original = points(vertices)?;
    let triangles = triangles(indices, original.len())?;
    let sources = points(pin_source)?;
    let destinations = points(pin_destination)?;
    anyhow::ensure!(
        sources.len() == destinations.len(),
        "pin arrays differ in length"
    );

    let adjacency = adjacency(original.len(), &triangles);
    let mut constraints = HashMap::<usize, Point>::new();
    for (&source, &destination) in sources.iter().zip(&destinations) {
        let vertex = nearest(&original, source);
        constraints.insert(vertex, original[vertex].add(destination.sub(source)));
    }
    anchor_unconstrained_components(&adjacency, &original, &mut constraints);

    let mut deformed = initial_guess.to_vec();
    for (&vertex, &target) in &constraints {
        deformed[vertex] = target;
    }

    for _ in 0..ITERATIONS {
        let rotations = local_step(&original, &deformed, &adjacency);
        deformed = global_step(&original, &adjacency, &rotations, &constraints, &deformed);
    }
    Ok(deformed)
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

fn adjacency(vertex_count: usize, triangles: &[[usize; 3]]) -> Vec<Vec<usize>> {
    let mut sets = vec![HashSet::new(); vertex_count];
    for triangle in triangles {
        for (a, b) in [
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ] {
            sets[a].insert(b);
            sets[b].insert(a);
        }
    }
    sets.into_iter()
        .map(|set| set.into_iter().collect())
        .collect()
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
    adjacency: &[Vec<usize>],
    original: &[Point],
    constraints: &mut HashMap<usize, Point>,
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
            for &next in &adjacency[vertex] {
                if !visited[next] {
                    visited[next] = true;
                    stack.push(next);
                }
            }
        }
        if !constrained {
            constraints.insert(component[0], original[component[0]]);
        }
    }
}


fn local_step(original: &[Point], deformed: &[Point], adjacency: &[Vec<usize>]) -> Vec<(f64, f64)> {
    (0..original.len())
        .map(|i| {
            let mut cosine = 0.0;
            let mut sine = 0.0;
            for &j in &adjacency[i] {
                let p = original[i].sub(original[j]);
                let q = deformed[i].sub(deformed[j]);
                cosine += p.x * q.x + p.y * q.y;
                sine += p.x * q.y - p.y * q.x;
            }
            let norm = cosine.hypot(sine);
            if norm > 1.0e-12 {
                (cosine / norm, sine / norm)
            } else {
                (1.0, 0.0)
            }
        })
        .collect()
}

fn rotate(point: Point, rotation: (f64, f64)) -> Point {
    Point {
        x: rotation.0 * point.x - rotation.1 * point.y,
        y: rotation.1 * point.x + rotation.0 * point.y,
    }
}

fn global_step(
    original: &[Point],
    adjacency: &[Vec<usize>],
    rotations: &[(f64, f64)],
    constraints: &HashMap<usize, Point>,
    previous: &[Point],
) -> Vec<Point> {
    let free = (0..original.len())
        .filter(|i| !constraints.contains_key(i))
        .collect::<Vec<_>>();
    let free_index = free
        .iter()
        .enumerate()
        .map(|(row, &vertex)| (vertex, row))
        .collect::<HashMap<_, _>>();
    if free.is_empty() {
        return previous.to_vec();
    }
    let mut bx = vec![0.0; free.len()];
    let mut by = vec![0.0; free.len()];
    for (row, &i) in free.iter().enumerate() {
        for &j in &adjacency[i] {
            let edge = original[i].sub(original[j]);
            let rhs = rotate(edge, rotations[i])
                .add(rotate(edge, rotations[j]))
                .scale(0.5);
            bx[row] += rhs.x;
            by[row] += rhs.y;
            if let Some(target) = constraints.get(&j) {
                bx[row] += target.x;
                by[row] += target.y;
            }
        }
    }
    let apply = |input: &[f64]| -> Vec<f64> {
        free.iter()
            .map(|&i| {
                let mut value = adjacency[i].len() as f64 * input[free_index[&i]];
                for &j in &adjacency[i] {
                    if let Some(&column) = free_index.get(&j) {
                        value -= input[column];
                    }
                }
                value
            })
            .collect()
    };
    let initial_x = free.iter().map(|&i| previous[i].x).collect::<Vec<_>>();
    let initial_y = free.iter().map(|&i| previous[i].y).collect::<Vec<_>>();
    let x = conjugate_gradient(&apply, &bx, initial_x);
    let y = conjugate_gradient(&apply, &by, initial_y);
    let mut result = previous.to_vec();
    for (row, &vertex) in free.iter().enumerate() {
        result[vertex] = Point {
            x: x[row],
            y: y[row],
        };
    }
    for (&vertex, &target) in constraints {
        result[vertex] = target;
    }
    result
}

fn conjugate_gradient(
    apply: &impl Fn(&[f64]) -> Vec<f64>,
    rhs: &[f64],
    mut x: Vec<f64>,
) -> Vec<f64> {
    let mut residual = subtract(rhs, &apply(&x));
    let mut direction = residual.clone();
    let mut residual_norm = dot(&residual, &residual);
    for _ in 0..CG_ITERATIONS.min(rhs.len().saturating_mul(4).max(1)) {
        if residual_norm.sqrt() < TOLERANCE {
            break;
        }
        let applied = apply(&direction);
        let denominator = dot(&direction, &applied);
        if denominator.abs() < 1.0e-20 {
            break;
        }
        let alpha = residual_norm / denominator;
        for i in 0..x.len() {
            x[i] += alpha * direction[i];
            residual[i] -= alpha * applied[i];
        }
        let next_norm = dot(&residual, &residual);
        let beta = next_norm / residual_norm;
        for i in 0..direction.len() {
            direction[i] = residual[i] + beta * direction[i];
        }
        residual_norm = next_norm;
    }
    x
}

fn subtract(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(a, b)| a - b).collect()
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

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
