use std::{result, usize};

use super::utils::round_to_n_decimal_places;
use nalgebra::{storage::Storage, Dim, Scalar};
use nalgebra::{
    ArrayStorage, Const, DMatrix, DimName, Matrix, Matrix2, SMatrix, SVector, Vector, Vector2,
    Vector4,
};
use rand::prelude::*;
use std::fmt;

pub fn s2d<const N: usize>(matrix: SMatrix<f64, N, N>) -> DMatrix<f64> {
    // This creates a DMatrix from the row slice of the SMatrix.
    DMatrix::from_row_slice(matrix.nrows(), matrix.ncols(), matrix.as_slice())
}

#[derive(Clone, Debug)]
pub struct Basis {
    bases_0: Vec<Vector2<f64>>,
    bases_1: Vec<Vector2<f64>>,
}

impl Basis {
    pub fn new(n: usize) -> Self {
        let mut bases_0: Vec<Vector2<f64>> = vec![];
        let mut bases_1: Vec<Vector2<f64>> = vec![];

        for _i in 0..n {
            bases_0.push(Vector2::new(1.0, 0.0));
            bases_1.push(Vector2::new(0.0, 1.0));
        }
        Self { bases_0, bases_1 }
    }
    pub fn get_basis(&self, index: usize) -> (Vector2<f64>, Vector2<f64>) {
        (self.bases_0[index], self.bases_1[index])
    }
    pub fn rotate_to_angle(&mut self, teta: f64, index: usize) -> Self {
        {
            // convert teta angle to radians
            // PI / 2.0 is 90 degres in radian
            let radians = teta.to_radians();

            let new_basis =
                Matrix2::new(radians.cos(), -radians.sin(), radians.sin(), radians.cos());

            self.bases_0[index] = new_basis.column(0).into();
            self.bases_1[index] = new_basis.column(1).into();

            self.clone()
        }
    }
}

#[derive(Clone)]
pub struct EntangledParticleN<const N: usize> {
    // v*w = r a_0*b_0 + s a_0*b_1 + t a_1*b_0 + u a_1*b_1
    pub basis: Basis,
    state: SVector<f64, N>,
}

impl<const N: usize> fmt::Debug for EntangledParticleN<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "EntangledParticleN {{")?;

        writeln!(f, "  state:")?;
        for (i, amp) in self.state.iter().enumerate() {
            let bits = format!("{:0width$b}", i, width = (N as f64).log2().ceil() as usize);
            writeln!(f, "    Amplitude {}: |{}⟩ has amplitude {}", i, bits, amp)?;
        }

        writeln!(f, "}}")
    }
}

impl<const N: usize> EntangledParticleN<N> {
    pub fn change_state_by_matrix(&mut self, matrix: SMatrix<f64, N, N>) {
        assert_eq!(N, self.state.len());
        let res = matrix * self.state;
        self.state = res;
        self.check_params();
    }

    pub fn get_no_particles(&self) -> usize {
        self.basis.bases_0.len()
    }

    pub fn get_params(&self) -> &SVector<f64, N> {
        self.check_params();
        &self.state
    }

    pub fn get_basis(&self) -> Basis {
        self.basis.clone()
    }
    pub fn new(state: SVector<f64, N>) -> Self {
        let res = Self {
            basis: Basis::new((usize::BITS - 1 - N.leading_zeros()) as usize),
            state,
        };
        res.check_params();
        res
    }
    pub fn check_params(&self) -> bool {
        let mut square_sum: f64 = 0.0;
        for i in self.state.iter() {
            square_sum += i.powi(2);
        }
        if round_to_n_decimal_places(square_sum, 5) != 1.0 {
            dbg!(self);
            panic!("params dont check out");
        }

        true
    }

    pub fn swap_basis(&mut self, angle: f64, index: usize) {
        let rotated_basis = self.get_basis().rotate_to_angle(angle, index);

        let alice_basis = rotated_basis.get_basis(index);

        self.basis.bases_0[index] = alice_basis.0;
        self.basis.bases_1[index] = alice_basis.1;
    }

    pub fn measure(&mut self, index: usize) -> bool {
        let n: usize = self.get_no_particles();

        let mut p_matrices: Vec<(Matrix2<f64>, Matrix2<f64>)> = vec![];
        // go through all particles
        for i in 0..n {
            // case when the particle is the one we wabt to measure
            if i == index {
                let prtcl_to_measure_0 = self.basis.bases_0[i] * self.basis.bases_0[i].transpose();
                let prtcl_to_measure_1 = self.basis.bases_1[i] * self.basis.bases_1[i].transpose();
                p_matrices.push((prtcl_to_measure_0, prtcl_to_measure_1));
            }
            // if its not the particle we want to measure we use identiy matrix
            else {
                p_matrices.push((Matrix2::identity(), Matrix2::identity()));
            }
        }

        let mut kronecker = s2d(p_matrices[0].0);

        for i in 1..self.get_no_particles() {
            kronecker = kronecker.kronecker(&p_matrices[i].0);
        }
        let mut new_state = kronecker * self.state;
        let probability = new_state.norm_squared();

        //sample
        let obs = rand::rng().random_bool(round_to_n_decimal_places(probability, 5));
        let mut norm = probability.sqrt();

        if !obs {
            norm = (1.0 - probability).sqrt();

            let mut kronecker_1 = s2d(p_matrices[0].1);

            for i in 1..self.get_no_particles() {
                kronecker_1 = kronecker_1.kronecker(&p_matrices[i].1);
            }

            new_state = kronecker_1 * self.state;
        }

        let post_state = new_state / norm;

        self.state = SMatrix::from_column_slice(post_state.as_slice());

        obs
    }
}

mod tests {

    use super::*;

    #[test]
    fn test_just_rotate() {
        let mut prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(
            0.5f64.sqrt(),
            0.0,
            0.0,
            0.5f64.sqrt(),
        ));
        dbg!(&prtcl.basis);
        prtcl.swap_basis(90.0, 0);
        dbg!(&prtcl.basis);
    }

    #[test]
    fn test_rotate_same_direction() {
        for i in 0..10 {
            let mut prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(
                0.5f64.sqrt(),
                0.0,
                0.0,
                0.5f64.sqrt(),
            ));

            prtcl.swap_basis(0.0, 0);
            prtcl.swap_basis(0.0, 1);

            let measurement_1 = prtcl.measure(0);
            let measurement_2 = prtcl.measure(1);
            assert_eq!(measurement_1, measurement_2);
        }

        for i in 0..10 {
            let mut prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(
                0.5f64.sqrt(),
                0.0,
                0.0,
                0.5f64.sqrt(),
            ));
            prtcl.swap_basis(0.0, 0);
            prtcl.swap_basis(0.0, 1);

            let measurement_1 = prtcl.measure(0);
            let measurement_2 = prtcl.measure(1);
            assert_eq!(measurement_1, measurement_2);
        }
    }
    #[test]
    fn test_rotate_not_same_direction() {
        let reps = 1000;
        let mut alice_count = 0;
        let mut bob_count = 0;
        let mut agreements = 0;
        for _i in 0..reps {
            let mut prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(
                0.5f64.sqrt(),
                0.0,
                0.0,
                0.5f64.sqrt(),
            ));

            prtcl.swap_basis(0.0, 0);
            prtcl.swap_basis(90.0 * (2.0 / 3.0), 1);

            let alice: bool = prtcl.measure(0);
            let bob: bool = prtcl.measure(1);

            if alice {
                alice_count += 1;
            };
            if bob {
                bob_count += 1;
            };

            if alice == bob {
                agreements += 1;
            }
        }

        dbg!(agreements);
        dbg!("alice %: {}", alice_count as f64 / reps as f64);
        dbg!("bob %: {}", bob_count as f64 / reps as f64);
        dbg!("agreements % {}", agreements as f64 / reps as f64);
        dbg!("expected 1/4");
    }
}
