//! Dependency-free SPD backend; geometry and boundary elimination stay in the caller.
use std::{
    cell::RefCell,
    cmp::Reverse,
    collections::{BTreeMap, BinaryHeap},
    rc::Rc,
};
#[derive(Clone, Copy)]
pub(crate) struct SolveOptions {
    pub max_iterations: usize,
    pub relative_tolerance: f64,
}
impl Default for SolveOptions {
    fn default() -> Self {
        Self {
            max_iterations: 96,
            relative_tolerance: 1e-8,
        }
    }
}
#[cfg(test)]
thread_local! {
    static FACTORIZATIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static ITERATIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}
#[cfg(test)]
pub(crate) fn take_work() -> [usize; 2] {
    [
        ITERATIONS.with(|n| n.replace(0)),
        FACTORIZATIONS.with(|n| n.replace(0)),
    ]
}
/// A_ii = diagonal[i], A_ij = -rows[i][j], using compact free indices.
/// A future sparse library can replace this type without touching deformation.
pub(crate) struct PreparedSystem {
    rows: Vec<Vec<(usize, f64)>>,
    diagonal: Vec<f64>,
    factor: Option<Rc<LinearFactor>>,
    preconditioner: Option<Rc<IncompleteCholesky>>,
    options: SolveOptions,
}
impl PreparedSystem {
    pub fn new(rows: Vec<Vec<(usize, f64)>>, diagonal: Vec<f64>, options: SolveOptions) -> Self {
        let (factor, preconditioner) = cached_factor(&rows, &diagonal);
        Self {
            rows,
            diagonal,
            factor,
            preconditioner,
            options,
        }
    }
    pub fn solve(&self, rhs: &[f64], initial: Vec<f64>) -> Vec<f64> {
        if let Some(factor) = &self.factor {
            factor.solve(rhs)
        } else {
            conjugate_gradient_steps(
                &|input| self.apply(input),
                rhs,
                initial,
                &self.diagonal,
                self.options.max_iterations,
                self.options.relative_tolerance,
                self.preconditioner.as_deref(),
            )
        }
    }
    fn apply(&self, input: &[f64]) -> Vec<f64> {
        self.rows
            .iter()
            .enumerate()
            .map(|(i, row)| {
                self.diagonal[i] * input[i] - row.iter().map(|&(j, w)| w * input[j]).sum::<f64>()
            })
            .collect()
    }
    #[cfg(test)]
    pub fn uses_sparse_factor(&self) -> bool {
        matches!(self.factor.as_deref(), Some(LinearFactor::Sparse(_)))
    }
}

// The factors depend on the rest mesh and constrained indices, never target
// positions. Keep a bounded per-thread cache so a drag only updates right sides.
struct FactorEntry {
    rows: Vec<Vec<(usize, f64)>>,
    diagonal: Vec<f64>,
    factor: Option<Rc<LinearFactor>>,
    preconditioner: Option<Rc<IncompleteCholesky>>,
}
thread_local! {static FACTORS: RefCell<Vec<FactorEntry>> = const { RefCell::new(Vec::new()) };}
fn cached_factor(
    rows: &[Vec<(usize, f64)>],
    diagonal: &[f64],
) -> (Option<Rc<LinearFactor>>, Option<Rc<IncompleteCholesky>>) {
    FACTORS.with(|cache| {
        let mut cache = cache.borrow_mut();
        if let Some(i) = cache
            .iter()
            .position(|e| e.diagonal == diagonal && e.rows == rows)
        {
            let entry = cache.remove(i);
            let result = (entry.factor.clone(), entry.preconditioner.clone());
            cache.push(entry);
            return result;
        }
        #[cfg(test)]
        FACTORIZATIONS.with(|n| n.set(n.get() + 1));
        let n = diagonal.len();
        let factor = if n > 0 && n <= 384 {
            let mut matrix = vec![0.; n * n];
            for i in 0..n {
                matrix[i * n + i] = diagonal[i];
                for &(j, w) in &rows[i] {
                    matrix[i * n + j] = -w;
                }
            }
            cholesky(matrix, n).map(LinearFactor::Dense)
        } else {
            SparseFactor::new(rows, diagonal).map(LinearFactor::Sparse)
        }
        .map(Rc::new);
        let preconditioner = if factor.is_none() {
            IncompleteCholesky::new(rows, diagonal).map(Rc::new)
        } else {
            None
        };
        if cache.len() == 6 {
            cache.remove(0);
        }
        cache.push(FactorEntry {
            rows: rows.to_vec(),
            diagonal: diagonal.to_vec(),
            factor: factor.clone(),
            preconditioner: preconditioner.clone(),
        });
        (factor, preconditioner)
    })
}
enum LinearFactor {
    Dense(Vec<f64>),
    Sparse(SparseFactor),
}
impl LinearFactor {
    fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        match self {
            Self::Dense(a) => cholesky_solve(a, rhs),
            Self::Sparse(a) => a.solve(rhs),
        }
    }
}
// Minimum-degree sparse Cholesky. Elimination and fill depend only on matrix
// sparsity; an extreme target does not request a more expensive solve.
struct SparseFactor {
    columns: Vec<(usize, f64, Vec<(usize, f64)>)>,
}
impl SparseFactor {
    fn new(rows: &[Vec<(usize, f64)>], diagonal: &[f64]) -> Option<Self> {
        let n = rows.len();
        let mut work = rows
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&(j, w)| (j, -w))
                    .collect::<BTreeMap<_, _>>()
            })
            .collect::<Vec<_>>();
        let mut d = diagonal.to_vec();
        let mut active = vec![true; n];
        let mut queue = (0..n)
            .map(|i| Reverse((work[i].len(), i)))
            .collect::<BinaryHeap<_>>();
        let mut columns = Vec::with_capacity(n);
        let mut entries = 0;
        while columns.len() < n {
            let Reverse((degree, i)) = queue.pop()?;
            if !active[i] || work[i].len() != degree {
                continue;
            }
            let pivot = d[i];
            if !pivot.is_finite() || pivot <= 0. {
                return None;
            }
            let root = pivot.sqrt();
            let neighbors = work[i]
                .iter()
                .map(|(&j, &v)| (j, v / root))
                .collect::<Vec<_>>();
            entries += neighbors.len();
            if entries > 2_000_000 {
                return None;
            }
            active[i] = false;
            for (k, &(j, l)) in neighbors.iter().enumerate() {
                d[j] -= l * l;
                for &(q, r) in &neighbors[k + 1..] {
                    let value = work[j].get(&q).copied().unwrap_or(0.) - l * r;
                    work[j].insert(q, value);
                    work[q].insert(j, value);
                }
                work[j].remove(&i);
            }
            work[i].clear();
            for &(j, _) in &neighbors {
                queue.push(Reverse((work[j].len(), j)));
            }
            columns.push((i, root, neighbors));
        }
        Some(Self { columns })
    }
    fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        let mut x = rhs.to_vec();
        for (i, d, column) in &self.columns {
            x[*i] /= d;
            for &(j, l) in column {
                x[j] -= l * x[*i];
            }
        }
        for (i, d, column) in self.columns.iter().rev() {
            for &(j, l) in column {
                x[*i] -= l * x[j];
            }
            x[*i] /= d;
        }
        x
    }
}

fn cholesky(mut a: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    for i in 0..n {
        for j in 0..=i {
            let mut value = a[i * n + j];
            for k in 0..j {
                value -= a[i * n + k] * a[j * n + k];
            }
            if i == j {
                if value <= 0.0 || !value.is_finite() {
                    return None;
                }
                a[i * n + j] = value.sqrt();
            } else {
                a[i * n + j] = value / a[j * n + j];
            }
        }
    }
    Some(a)
}
fn cholesky_solve(a: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    let mut x = rhs.to_vec();
    for i in 0..n {
        let mut value = x[i];
        for j in 0..i {
            value -= a[i * n + j] * x[j];
        }
        x[i] = value / a[i * n + i];
    }
    for i in (0..n).rev() {
        let mut value = x[i];
        for j in i + 1..n {
            value -= a[j * n + i] * x[j];
        }
        x[i] = value / a[i * n + i];
    }
    x
}

// IC(0) keeps the sparse matrix pattern. A diagonal shift is used only in
// the preconditioner when an incomplete factor would have a nonpositive pivot;
// the actual system and its residual remain unchanged.
struct IncompleteCholesky {
    lower: Vec<Vec<(usize, f64)>>,
    diagonal: Vec<f64>,
}
impl IncompleteCholesky {
    fn new(rows: &[Vec<(usize, f64)>], diagonal: &[f64]) -> Option<Self> {
        for shift in [0., 0.001, 0.01, 0.1, 1., 10.] {
            let mut lower = Vec::<Vec<(usize, f64)>>::new();
            let mut pivots = Vec::<f64>::new();
            let mut valid = true;
            for (i, input) in rows.iter().enumerate() {
                let mut row = input
                    .iter()
                    .filter(|(j, _)| *j < i)
                    .map(|&(j, w)| (j, -w))
                    .collect::<Vec<_>>();
                row.sort_unstable_by_key(|&(j, _)| j);
                for k in 0..row.len() {
                    let j = row[k].0;
                    let (mut a, mut b, mut product) = (0, 0, 0.);
                    while a < k && b < lower[j].len() {
                        let (x, v) = row[a];
                        let (y, w) = lower[j][b];
                        if x == y {
                            product += v * w;
                            a += 1;
                            b += 1;
                        } else if x < y {
                            a += 1;
                        } else {
                            b += 1;
                        }
                    }
                    row[k].1 = (row[k].1 - product) / pivots[j];
                }
                let pivot =
                    diagonal[i] * (1. + shift) - row.iter().map(|&(_, v)| v * v).sum::<f64>();
                if !pivot.is_finite() || pivot <= diagonal[i] * 1e-10 {
                    valid = false;
                    break;
                }
                lower.push(row);
                pivots.push(pivot.sqrt());
            }
            if valid {
                return Some(Self {
                    lower,
                    diagonal: pivots,
                });
            }
        }
        None
    }
    fn apply(&self, rhs: &[f64]) -> Vec<f64> {
        let mut result = rhs.to_vec();
        for i in 0..result.len() {
            result[i] = (result[i]
                - self.lower[i]
                    .iter()
                    .map(|&(j, l)| l * result[j])
                    .sum::<f64>())
                / self.diagonal[i];
        }
        for i in (0..result.len()).rev() {
            result[i] /= self.diagonal[i];
            for &(j, l) in &self.lower[i] {
                result[j] -= l * result[i];
            }
        }
        result
    }
}

fn conjugate_gradient_steps(
    apply: &impl Fn(&[f64]) -> Vec<f64>,
    rhs: &[f64],
    mut x: Vec<f64>,
    diagonal: &[f64],
    limit: usize,
    tolerance: f64,
    preconditioner: Option<&IncompleteCholesky>,
) -> Vec<f64> {
    let mut residual = subtract(rhs, &apply(&x));
    let precondition = |r: &[f64]| {
        if let Some(factor) = preconditioner {
            factor.apply(r)
        } else {
            r.iter().zip(diagonal).map(|(r, d)| r / d).collect()
        }
    };
    let mut z = precondition(&residual);
    let mut direction = z.clone();
    let mut residual_norm = dot(&residual, &z);
    let tolerance = tolerance * tolerance * dot(rhs, rhs).max(1.0);
    for _ in 0..limit.min(rhs.len().saturating_mul(4).max(1)) {
        if dot(&residual, &residual) < tolerance {
            break;
        }
        let applied = apply(&direction);
        let denominator = dot(&direction, &applied);
        if denominator.abs() < 1.0e-20 {
            break;
        }
        #[cfg(test)]
        ITERATIONS.with(|n| n.set(n.get() + 1));
        let alpha = residual_norm / denominator;
        for i in 0..x.len() {
            x[i] += alpha * direction[i];
            residual[i] -= alpha * applied[i];
        }
        z = precondition(&residual);
        let next_norm = dot(&residual, &z);
        let beta = next_norm / residual_norm;
        for i in 0..direction.len() {
            direction[i] = z[i] + beta * direction[i];
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
    fn iterative_backend_matches_direct_solution_and_respects_its_budget() {
        let rows = vec![
            vec![(1, 1.), (3, 0.7)],
            vec![(0, 1.), (2, 2.)],
            vec![(1, 2.), (3, 0.5)],
            vec![(0, 0.7), (2, 0.5)],
        ];
        let diagonal = vec![3., 4., 5., 6.];
        let mut system = PreparedSystem::new(
            rows,
            diagonal,
            SolveOptions {
                max_iterations: 32,
                relative_tolerance: 1e-12,
            },
        );
        let exact = [1.2, -2., 3.4, 0.7];
        let rhs = system.apply(&exact);
        let direct = system.solve(&rhs, vec![0.; 4]);
        system.factor = None;
        system.preconditioner = None;
        let iterative = system.solve(&rhs, vec![0.; 4]);
        for ((a, b), c) in direct.iter().zip(&iterative).zip(exact) {
            assert!((a - c).abs() < 1e-11 && (b - c).abs() < 1e-11);
        }
        system.options.max_iterations = 0;
        assert_eq!(system.solve(&rhs, vec![0.; 4]), vec![0.; 4]);
    }
    #[test]
    fn changed_right_hand_sides_reuse_the_same_factor() {
        let rows = vec![vec![(1, 0.31)], vec![(0, 0.31)]];
        let diagonal = vec![1.123, 2.456];
        let a = PreparedSystem::new(rows.clone(), diagonal.clone(), SolveOptions::default());
        take_work();
        let b = PreparedSystem::new(rows, diagonal, SolveOptions::default());
        for rhs in [vec![2., 3.], vec![-10., 27.]] {
            let x = b.solve(&rhs, vec![0.; 2]);
            for (actual, expected) in a.apply(&x).iter().zip(rhs) {
                assert!((actual - expected).abs() < 1e-12);
            }
        }
        assert_eq!(take_work(), [0, 0]);
    }
}
