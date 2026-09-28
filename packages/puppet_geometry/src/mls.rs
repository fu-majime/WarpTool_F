//! Point-based rigid MLS (Schaefer, McPhail, Warren, 2006).
//! Explicit user scale is prescribed, never estimated from stretched targets.
use crate::{math::Point, pose::Pose};

struct Control {
    source: Point,
    target: Point,
    scale: f64,
}

pub(crate) struct RigidMls {
    controls: Vec<Control>,
    exponent: f64,
    identity: bool,
}

impl RigidMls {
    pub fn new(sources: &[f64], targets: &[f64], poses: &[Pose], exponent: f64) -> Self {
        let mut poses = poses.to_vec();
        if sources.is_empty() && !poses.is_empty() && poses.iter().all(|p| !p.moving) {
            place_free_frames(&mut poses);
        }
        // Real pins precede frame samples, so coincident observations cannot
        // override an exact positional target.
        let mut controls = sources
            .chunks_exact(2)
            .zip(targets.chunks_exact(2))
            .map(|(p, q)| Control {
                source: Point { x: p[0], y: p[1] },
                target: Point { x: q[0], y: q[1] },
                scale: poses
                    .iter()
                    .find(|h| h.source == (p[0], p[1]))
                    .map_or(1., |h| h.scale),
            })
            .collect::<Vec<_>>();
        for p in poses.iter().filter(|p| p.moving || p.rotate || p.resize) {
            let angle = p.natural_angle + if p.rotate { p.angle } else { 0. };
            let (s, c) = angle.sin_cos();
            // A symmetric frame encodes rotation/scale as ordinary MLS point
            // observations. No mesh territories, silhouette sectors or post-warp.
            for (x, y) in [
                (p.radius, 0.),
                (-p.radius, 0.),
                (0., p.radius),
                (0., -p.radius),
            ] {
                controls.push(Control {
                    source: Point {
                        x: p.source.0 + x,
                        y: p.source.1 + y,
                    },
                    target: Point {
                        x: p.center.0 + p.scale * (c * x - s * y),
                        y: p.center.1 + p.scale * (s * x + c * y),
                    },
                    scale: p.scale,
                });
            }
        }
        let identity = controls
            .iter()
            .all(|c| c.source == c.target && c.scale == 1.);
        Self {
            controls,
            exponent,
            identity,
        }
    }

    pub fn evaluate(&self, point: Point) -> Point {
        if self.identity {
            return point;
        }
        let distances = self
            .controls
            .iter()
            .map(|c| {
                let d = point.sub(c.source);
                d.x.hypot(d.y)
            })
            .collect::<Vec<_>>();
        if let Some(i) = distances.iter().position(|&d| d == 0.) {
            return self.controls[i].target;
        }
        let nearest = distances.iter().copied().fold(f64::INFINITY, f64::min);
        // Normalize before exponentiation to avoid overflow near a pin.
        let weights = distances
            .iter()
            .map(|&d| (nearest / d).powf(2. * self.exponent))
            .collect::<Vec<_>>();
        let sum = weights.iter().sum::<f64>();
        let (mut p, mut q, mut scale) = (Point::default(), Point::default(), 0.);
        for (c, &w) in self.controls.iter().zip(&weights) {
            let w = w / sum;
            p = p.add(c.source.scale(w));
            q = q.add(c.target.scale(w));
            scale += w * c.scale;
        }
        let (mut a, mut b) = (0., 0.);
        for (c, &w) in self.controls.iter().zip(&weights) {
            let u = c.source.sub(p);
            let v = c.target.sub(q);
            a += w * (u.x * v.x + u.y * v.y);
            b += w * (u.x * v.y - u.y * v.x);
        }
        let norm = a.hypot(b);
        let (c, s) = if norm > 0. {
            (a / norm, b / norm)
        } else {
            (1., 0.)
        };
        let d = point.sub(p);
        q.add(
            Point {
                x: c * d.x - s * d.y,
                y: s * d.x + c * d.y,
            }
            .scale(scale),
        )
    }
}

/// Least-squares integration of frame vectors on the complete control graph:
/// min sum_ij |q_i-q_j - (F_i+F_j)(p_i-p_j)/2|^2, with q_0 fixed.
/// The complete-graph Laplacian has a closed-form inverse modulo translation.
/// This uses all frames symmetrically, without a selected tree or angle unwrap.
fn place_free_frames(poses: &mut [Pose]) {
    let frame = |p: &Pose| {
        let angle = p.natural_angle + if p.rotate { p.angle } else { 0. };
        Point {
            x: p.scale * angle.cos(),
            y: p.scale * angle.sin(),
        }
    };
    let source = |p: &Pose| Point {
        x: p.source.0,
        y: p.source.1,
    };
    let product = |z: Point, p: Point| Point {
        x: z.x * p.x - z.y * p.y,
        y: z.y * p.x + z.x * p.y,
    };
    let mean = |f: fn(&Pose) -> Point| {
        poses
            .iter()
            .map(f)
            .fold(Point::default(), Point::add)
            .scale(1. / poses.len() as f64)
    };
    let center = mean(source);
    let rotation = mean(frame);
    let base = source(&poses[0]);
    let reference = product(frame(&poses[0]), base.sub(center));
    let anchor = Point {
        x: poses[0].center.0,
        y: poses[0].center.1,
    };
    for p in poses {
        let offset = product(frame(p), source(p).sub(center))
            .sub(reference)
            .add(product(rotation, source(p).sub(base)))
            .scale(0.5);
        let q = anchor.add(offset);
        p.center = (q.x, q.y);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn free_frame_centers_satisfy_the_complete_graph_normal_equations() {
        let mut poses = crate::pose::parse(&[
            -2., 0., 3., 4., 0., 1., 0.2, 1., 2., 0., 1., 1., 1., 2., 1., 2., 0., 1., -0.7, 1.2,
            2., 0., 1., 1., 3., -1., 3., -1., 0., 1., 1.1, 0.8, 2., 0., 1., 1.,
        ])
        .unwrap();
        place_free_frames(&mut poses);
        assert_eq!(poses[0].center, (3., 4.));
        for p in &poses[1..] {
            let mut residual = Point::default();
            for q in &poses {
                let x = p.source.0 - q.source.0;
                let y = p.source.1 - q.source.1;
                let rotate = |h: &Pose| Point {
                    x: h.scale * (h.angle.cos() * x - h.angle.sin() * y),
                    y: h.scale * (h.angle.sin() * x + h.angle.cos() * y),
                };
                let desired = rotate(p).add(rotate(q)).scale(0.5);
                residual = residual.add(
                    Point {
                        x: p.center.0 - q.center.0,
                        y: p.center.1 - q.center.1,
                    }
                    .sub(desired),
                );
            }
            assert!(residual.x.hypot(residual.y) < 1e-12);
        }
    }

    #[test]
    fn nonrigid_targets_still_produce_the_rigid_least_squares_fit() {
        // All three observations are equidistant from (0,0). Their centered
        // covariance has trace 11/3 and skew 1/3; SO(2), not similarity/affine.
        let field = RigidMls::new(
            &[-1., 0., 1., 0., 0., 1.],
            &[-1., 0., 2., 0., 0., 1.],
            &[],
            1.,
        );
        let q = field.evaluate(Point { x: 0., y: 0. });
        let norm = 122_f64.sqrt();
        assert!((q.x - (1. / 3. + 1. / (3. * norm))).abs() < 1e-12);
        assert!((q.y - (1. / 3. - 11. / (3. * norm))).abs() < 1e-12);
        let stretched = RigidMls::new(&[-1., 0., 1., 0.], &[-2., 0., 2., 0.], &[], 1.);
        assert_eq!(
            stretched.evaluate(Point { x: 0., y: 1. }),
            Point { x: 0., y: 1. }
        );
    }

    #[test]
    fn rigid_motion_and_real_pin_interpolation_are_exact_to_roundoff() {
        let source = [-2., -1., 3., 2., 0., 4.];
        for angle in [0., 0.7, std::f64::consts::PI, 7.] {
            let transform = |p: Point| Point {
                x: angle.cos() * p.x - angle.sin() * p.y + 5.,
                y: angle.sin() * p.x + angle.cos() * p.y - 2.,
            };
            let target = source
                .chunks_exact(2)
                .flat_map(|p| {
                    let q = transform(Point { x: p[0], y: p[1] });
                    [q.x, q.y]
                })
                .collect::<Vec<_>>();
            let field = RigidMls::new(&source, &target, &[], 1.);
            for p in [
                Point { x: -3., y: 8. },
                Point { x: 0., y: 0. },
                Point { x: 0.2, y: -1.5 },
            ] {
                let q = field.evaluate(p);
                let e = transform(p);
                assert!((q.x - e.x).hypot(q.y - e.y) < 1e-12);
            }
            for (p, q) in source.chunks_exact(2).zip(target.chunks_exact(2)) {
                assert_eq!(
                    field.evaluate(Point { x: p[0], y: p[1] }),
                    Point { x: q[0], y: q[1] }
                );
            }
        }
    }
}
