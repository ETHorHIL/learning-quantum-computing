use crate::module::entangled_particle_n::EntangledParticleN;
use crate::module::gates::basic::*;
use crate::module::{
    gates::{
        self,
        basic::{self, SingleInputGate},
    },
    utils::round_to_n_decimal_places,
};
use nalgebra::SVector;
use nalgebra::{ArrayStorage, Const};

mod tests {

    use std::os::macos::raw::stat;

    use nalgebra::{Matrix, Matrix2, Matrix3, Matrix4};

    use crate::module::entangled_particle_n;

    use super::*;

    fn deutsch(constant: bool) {
        let mut inputs = [0.0; 16].to_vec();
        inputs[1] = 1.0;

        let mut state = EntangledParticleN::new(SVector::<f64, 8>::from_vec(inputs));
        dbg!(&state);

        let hadamard2x2 = SingleInputGate::Hadamard.get_matrix::<2>(0);
        let hadamard8x8 = hadamard2x2.kronecker(&hadamard2x2).kronecker(&hadamard2x2);
        let p_0 = Matrix2::new(1.0, 0.0, 0.0, 0.0);
        let p_1 = Matrix2::new(0.0, 0.0, 0.0, 1.0);
        let identity_2x2 = Matrix2::identity();
        //let identity_4x4 = Matrix4::identity();
        let x = SingleInputGate::X.get_matrix::<2>(0);

        let f_balanced = p_0.kronecker(&identity_2x2).kronecker(&identity_2x2)
            + p_1.kronecker(&identity_2x2).kronecker(&x);

        let f_constant = Matrix::<f64, Const<8>, Const<8>, ArrayStorage<f64, 8, 8>>::identity();

        // hadamard to input
        // puts inputs into superposition
        state.change_state_by_matrix(hadamard8x8);

        // f constant will put the superposition state into state where all phases are equal. I.e. all amplitudes are either -1 or 1
        // f balanced will put the superposition into a state where amplitudes are alternatively -1 or 1 with the same amount of each occuring
        if constant {
            state.change_state_by_matrix(f_constant);
        } else {
            state.change_state_by_matrix(f_balanced);
        }

        // hadamard on input
        // f_constant: if the input is in an equal weighted superposition then applying hadamard will put the state back to pure state 00. So |00> will have amplitude 1.
        // f_balance: if the input is balanced then the hadamard gate will not put it back into a pure state. 00 will not have probability 1.
        state.change_state_by_matrix(hadamard8x8);

        state.measure(0);
        state.measure(1);
        state.measure(2);
        state.measure(3);
        state.measure(4);

        // if there is a non zero probability (amplitude =!0) of measuring |00> then the function must
        // be balanced

        let res = state.get_params().column(0)[1];
        if constant {
            assert_eq!(res, 1.0);
        } else {
            assert_ne!(res, 1.0);
        }
    }
    #[test]
    fn test_deutsch() {
        deutsch(true);
        deutsch(false);
    }
}
