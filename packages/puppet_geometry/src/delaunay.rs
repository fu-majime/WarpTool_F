//! Robust CDT; classify against conservative polygons, not raster centroids.
use aviutl2::anyhow::{self, Context as _};
use spade::{ConstrainedDelaunayTriangulation, HasPosition, Point2, Triangulation};
#[derive(Debug)]
pub struct ConstraintConflict(pub Vec<usize>);
impl std::fmt::Display for ConstraintConflict {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "intersecting contour constraints")
    }
}
impl std::error::Error for ConstraintConflict {}
#[derive(Clone, Copy)]
struct Vertex {
    position: Point2<f64>,
    index: usize,
}
impl HasPosition for Vertex {
    type Scalar = f64;
    fn position(&self) -> Point2<f64> {
        self.position
    }
}
/// Validate the complete boundary before sampling or triangulating the interior.
/// The original pixel-cell rings are the topology oracle. A conflicting
/// approximation is replaced once, without an arbitrary retry/tolerance budget.
pub fn validated_contours(
    raw: &[Vec<(f32, f32)>],
    mut candidate: Vec<Vec<(f32, f32)>>,
) -> anyhow::Result<Vec<Vec<(f32, f32)>>> {
    anyhow::ensure!(raw.len() == candidate.len(), "contour count mismatch");
    let mut restored = vec![false; raw.len()];
    loop {
        let mut points = Vec::new();
        let mut edges = Vec::new();
        let mut owners = Vec::new();
        for (id, ring) in candidate.iter().enumerate() {
            anyhow::ensure!(ring.len() >= 3, "degenerate contour");
            let base = points.len();
            points.extend(ring.iter().map(|&(x, y)| (x as f64, y as f64)));
            owners.extend(std::iter::repeat_n(id, ring.len()));
            edges.extend((0..ring.len()).map(|i| (base + i, base + (i + 1) % ring.len())));
        }
        match constrained_mesh(&points, &edges) {
            Ok(_) => return Ok(candidate),
            Err(error) => {
                let Some(conflict) = error.downcast_ref::<ConstraintConflict>() else {
                    return Err(error);
                };
                let mut changed = false;
                for &vertex in &conflict.0 {
                    let id = owners[vertex];
                    if !restored[id] {
                        // Zero tolerance removes only redundant collinear raster
                        // corners. It preserves islands, holes and cell coverage.
                        candidate[id] = crate::simplify::simplify(&raw[id], 0., f32::MAX, 0.);
                        restored[id] = true;
                        changed = true;
                    }
                }
                anyhow::ensure!(changed, "invalid source contour topology: {error}");
            }
        }
    }
}

pub fn triangulate(
    points: &[(f32, f32)],
    edges: &[(usize, usize)],
) -> anyhow::Result<Vec<[usize; 3]>> {
    triangulate_f64(
        &points
            .iter()
            .map(|&(x, y)| (x as f64, y as f64))
            .collect::<Vec<_>>(),
        edges,
    )
}
pub fn triangulate_f64(
    points: &[(f64, f64)],
    edges: &[(usize, usize)],
) -> anyhow::Result<Vec<[usize; 3]>> {
    let cdt = constrained_mesh(points, edges)?;
    Ok(cdt
        .inner_faces()
        .filter_map(|face| {
            let vertices = face.vertices();
            let p = face.positions();
            let x = (p[0].x + p[1].x + p[2].x) / 3.0;
            let y = (p[0].y + p[1].y + p[2].y) / 3.0;
            let mut inside = edges.is_empty();
            for &(a, b) in edges {
                let (ax, ay) = points[a];
                let (bx, by) = points[b];
                if (ay > y) != (by > y) && x < ax + (y - ay) * (bx - ax) / (by - ay) {
                    inside = !inside;
                }
            }
            inside.then(|| vertices.map(|v| v.data().index))
        })
        .collect())
}

fn constrained_mesh(
    points: &[(f64, f64)],
    edges: &[(usize, usize)],
) -> anyhow::Result<ConstrainedDelaunayTriangulation<Vertex>> {
    anyhow::ensure!(
        edges
            .iter()
            .all(|&(a, b)| a < points.len() && b < points.len()),
        "constraint index out of bounds"
    );
    let mut cdt = ConstrainedDelaunayTriangulation::<Vertex>::new();
    let handles = points
        .iter()
        .enumerate()
        .map(|(index, &(x, y))| {
            cdt.insert(Vertex {
                position: Point2::new(x, y),
                index,
            })
            .context("invalid mesh vertex")
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    for &(a, b) in edges {
        if handles[a] == handles[b] {
            continue;
        }
        if cdt.try_add_constraint(handles[a], handles[b]).is_empty() {
            let mut affected = vec![a, b];
            for edge in cdt.get_conflicting_edges_between_vertices(handles[a], handles[b]) {
                affected.extend(edge.vertices().map(|v| v.data().index));
            }
            return Err(ConstraintConflict(affected).into());
        }
    }
    Ok(cdt)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn crossed_approximation_is_repaired_before_interior_sampling() {
        let raw = vec![vec![(0., 0.), (4., 0.), (4., 4.), (0., 4.)]];
        let invalid = vec![vec![(0., 0.), (4., 4.), (4., 0.), (0., 4.)]];
        let valid = validated_contours(&raw, invalid).unwrap();
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];
        let triangles = triangulate(&valid[0], &edges).unwrap();
        let area = triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| valid[0][i]);
                ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)).abs() * 0.5
            })
            .sum::<f32>();
        assert_eq!(area, 16.);
    }
    #[test]
    fn invalid_source_and_indices_are_reported_without_a_retry_budget() {
        let bow = vec![vec![(0., 0.), (4., 4.), (4., 0.), (0., 4.)]];
        assert!(validated_contours(&bow, bow.clone()).is_err());
        assert!(triangulate(&[(0., 0.), (1., 0.), (0., 1.)], &[(0, 3)]).is_err());
    }
}
