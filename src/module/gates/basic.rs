use std::convert::identity;

use crate::module::entangled_particle_n::{s2d, EntangledParticleN};
use nalgebra::{ArrayStorage, Const, Matrix, Matrix2, Matrix4, SMatrix, SVector, Vector2, Vector4};

/*pub fn reverse_bellcircuit(mut prctl: EntangledParticle2) -> EntangledParticle2 {
    prctl.change_state_by_matrix(Gate::CNot_control_0.get_matrix());
    prctl.change_state_by_matrix(Gate::Hadamard_0.get_matrix());
    prctl
}*/

pub fn log(N: usize) -> usize {
    assert!(N.is_power_of_two(), "Input must be a power of 2");
    N.trailing_zeros() as usize
}
pub enum SingleInputGate {
    X,
    Z,
    Y,
    Hadamard,
}

pub enum TwoInputGate {
    CNot,
}

impl SingleInputGate {
    pub fn get_matrix<const N: usize>(&self, input_index: usize) -> SMatrix<f64, N, N> {
        let n: usize = log(N);
        let matrix_1 = Matrix2::new(1.0, 0.0, 0.0, 1.0);

        let gate_matrix: Matrix2<f64> = match self {
            SingleInputGate::X => Matrix2::new(0.0, 1.0, 1.0, 0.0),
            SingleInputGate::Z => Matrix2::new(1.0, 0.0, 0.0, -1.0),
            SingleInputGate::Y => Matrix2::new(0.0, 1.0, -1.0, 0.0),
            SingleInputGate::Hadamard => 1.0 / 2.0f64.sqrt() * Matrix2::new(1.0, 1.0, 1.0, -1.0),
        };

        let mut kronecker = if input_index == 0 {
            s2d(gate_matrix)
        } else {
            s2d(matrix_1)
        };

        for i in 1..n {
            if i == input_index {
                kronecker = kronecker.kronecker(&gate_matrix);
            } else {
                kronecker = kronecker.kronecker(&matrix_1);
            }
        }
        //dbg!(&kronecker);
        assert_eq!(kronecker.len(), N.pow(2));
        SMatrix::from_column_slice(kronecker.as_slice())
    }
}

impl TwoInputGate {
    pub fn get_matrix<const N: usize>(
        &self,
        control_index: usize,
        input_index: usize,
    ) -> SMatrix<f64, N, N> {
        let n = log(N);
        let identity = SMatrix::<f64, 2, 2>::identity();

        let gate_matrix: Matrix4<f64> = match self {
            TwoInputGate::CNot => {
                // inputs need to be adjacent for this code to work
                assert!(control_index.abs_diff(input_index) == 1);

                // Define the 2x2 projectors and gates:
                let proj0 = SMatrix::<f64, 2, 2>::new(1.0, 0.0, 0.0, 0.0);
                let proj1 = SMatrix::<f64, 2, 2>::new(0.0, 0.0, 0.0, 1.0);

                let x_gate = SMatrix::<f64, 2, 2>::new(0.0, 1.0, 1.0, 0.0);

                // Construct the 2-qubit CNOT operator:
                if control_index < input_index {
                    proj0.kronecker(&identity) + proj1.kronecker(&x_gate)
                } else {
                    proj1.kronecker(&identity) + proj0.kronecker(&x_gate)
                }
            }
        };

        let mut kronecker = if input_index == 0 || control_index == 0 {
            s2d(gate_matrix)
        } else {
            s2d(identity)
        };

        for i in 2..n {
            if i == input_index {
                kronecker = kronecker.kronecker(&gate_matrix);
            } else {
                kronecker = kronecker.kronecker(&identity);
            }
        }
        //dbg!(&kronecker);
        assert_eq!(kronecker.len(), N.pow(2));
        SMatrix::from_column_slice(kronecker.as_slice())
    }
}

#[cfg(test)]
mod tests {

    use crate::module::utils::round_to_n_decimal_places;

    use super::*;

    #[test]
    fn test_gates() {
        let mut prtcl = EntangledParticleN::new(SVector::<f64, 2>::new(0.0, 1.0));
        let hadamard = SingleInputGate::Hadamard.get_matrix::<2>(0);
        prtcl.change_state_by_matrix(hadamard);
        let params = prtcl.get_params();
        assert_eq!(
            [
                round_to_n_decimal_places(params[0], 5),
                round_to_n_decimal_places(params[1], 5)
            ],
            [
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
                round_to_n_decimal_places(-(0.5f64.sqrt()), 5)
            ]
        );

        let mut prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(1.0, 0.0, 0.0, 0.0));
        let hadamard = SingleInputGate::Hadamard.get_matrix::<4>(0);

        prtcl.change_state_by_matrix(hadamard);
        let params = prtcl.get_params();
        assert_eq!(
            [
                round_to_n_decimal_places(params[0], 5),
                round_to_n_decimal_places(params[1], 5),
                round_to_n_decimal_places(params[2], 5),
                round_to_n_decimal_places(params[3], 5)
            ],
            [
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
                round_to_n_decimal_places(0.0, 5),
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
                round_to_n_decimal_places(0.0, 5),
            ]
        );
        prtcl.change_state_by_matrix(hadamard);
        let params = prtcl.get_params();
        assert_eq!(
            [
                round_to_n_decimal_places(params[0], 5),
                round_to_n_decimal_places(params[1], 5),
                round_to_n_decimal_places(params[2], 5),
                round_to_n_decimal_places(params[3], 5)
            ],
            [
                round_to_n_decimal_places(1.0, 5),
                round_to_n_decimal_places(0.0, 5),
                round_to_n_decimal_places(0.0, 5),
                round_to_n_decimal_places(0.0, 5),
            ]
        );

        let mut prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(
            0.0,
            0.5f64.sqrt(),
            0.0,
            0.5f64.sqrt(),
        ));
        let cnot = TwoInputGate::CNot.get_matrix::<4>(0, 1);

        prtcl.change_state_by_matrix(cnot);
        let params = prtcl.get_params();
        assert_eq!(
            [
                round_to_n_decimal_places(params[0], 5),
                round_to_n_decimal_places(params[1], 5),
                round_to_n_decimal_places(params[2], 5),
                round_to_n_decimal_places(params[3], 5)
            ],
            [
                round_to_n_decimal_places(0.0, 5),
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
                round_to_n_decimal_places(0.0, 5),
            ]
        );
    }
    #[test]
    fn test_bell_circuit() {
        let mut prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(1.0, 0.0, 0.0, 0.0));
        //EntangledParticle2::new(1.0 / 2.0f64.sqrt(), 0.0, 0.0, 1.0 / 2.0f64.sqrt());
        let hadamard = SingleInputGate::Hadamard.get_matrix::<4>(0);
        let cnot = TwoInputGate::CNot.get_matrix::<4>(0, 1);
        dbg!(&hadamard);
        dbg!(&cnot);

        prtcl.change_state_by_matrix(hadamard);
        prtcl.change_state_by_matrix(cnot);

        let params = prtcl.get_params();
        assert_eq!(
            [
                round_to_n_decimal_places(params[0], 5),
                round_to_n_decimal_places(params[1], 5),
                round_to_n_decimal_places(params[2], 5),
                round_to_n_decimal_places(params[3], 5)
            ],
            [
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
                round_to_n_decimal_places(0.0, 5),
                round_to_n_decimal_places(0.0, 5),
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
            ]
        );

        let mut prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(0.0, 1.0, 0.0, 0.0));
        //EntangledParticle2::new(1.0 / 2.0f64.sqrt(), 0.0, 0.0, 1.0 / 2.0f64.sqrt());
        prtcl.change_state_by_matrix(hadamard);
        prtcl.change_state_by_matrix(cnot);

        let params = prtcl.get_params();
        assert_eq!(
            [
                round_to_n_decimal_places(params[0], 5),
                round_to_n_decimal_places(params[1], 5),
                round_to_n_decimal_places(params[2], 5),
                round_to_n_decimal_places(params[3], 5)
            ],
            [
                round_to_n_decimal_places(0.0, 5),
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
                round_to_n_decimal_places(0.0, 5),
            ]
        );

        let mut prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(0.0, 0.0, 1.0, 0.0));
        //EntangledParticle2::new(1.0 / 2.0f64.sqrt(), 0.0, 0.0, 1.0 / 2.0f64.sqrt());
        prtcl.change_state_by_matrix(hadamard);
        prtcl.change_state_by_matrix(cnot);

        let params = prtcl.get_params();
        assert_eq!(
            [
                round_to_n_decimal_places(params[0], 5),
                round_to_n_decimal_places(params[1], 5),
                round_to_n_decimal_places(params[2], 5),
                round_to_n_decimal_places(params[3], 5)
            ],
            [
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
                round_to_n_decimal_places(0.0, 5),
                round_to_n_decimal_places(0.0, 5),
                -1.0 * round_to_n_decimal_places(0.5f64.sqrt(), 5)
            ]
        );

        let mut prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(0.0, 0.0, 0.0, 1.0));
        //EntangledParticle2::new(1.0 / 2.0f64.sqrt(), 0.0, 0.0, 1.0 / 2.0f64.sqrt());
        prtcl.change_state_by_matrix(hadamard);
        prtcl.change_state_by_matrix(cnot);

        let params = prtcl.get_params();
        assert_eq!(
            [
                round_to_n_decimal_places(params[0], 5),
                round_to_n_decimal_places(params[1], 5),
                round_to_n_decimal_places(params[2], 5),
                round_to_n_decimal_places(params[3], 5)
            ],
            [
                round_to_n_decimal_places(0.0, 5),
                round_to_n_decimal_places(0.5f64.sqrt(), 5),
                -1.0 * round_to_n_decimal_places(0.5f64.sqrt(), 5),
                round_to_n_decimal_places(0.0, 5),
            ]
        );
    }
}
