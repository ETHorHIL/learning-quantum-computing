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

    use nalgebra::{DMatrix, Matrix, Matrix2, Matrix4, SMatrix};

    use crate::module::entangled_particle_n;

    use super::*;

    fn oracle() -> SMatrix<f64, 8, 8> {
        /*

        the oracle represents the "data" to be searched we have the following function
        f(0,0) = 0
        f(0,1) = 0
        f(1,0) = 1
        f(1,1) = 0

        we also have a controll qbit. So our if our basis is 000 ... 111 then the oracle can be represented as
        a matrix that is unity except for column 5 and 6 are flipped

        for a different encoding, e.g. the first entry being zero like this:
        f(0,0) = 1
        f(0,1) = 0
        f(1,0) = 0
        f(1,1) = 0
        swap column 1 and 2.
        For second bit swap column 3 and 4. Third bit 5 and six (as above), fourth bit 7 and 8.
        */
        let mut m = DMatrix::<f64>::identity(8, 8);
        m.swap_columns(4, 5);

        dbg!("oracle:");
        println!("{}", m);
        SMatrix::from_column_slice(m.as_slice())
    }

    fn amplifier() -> SMatrix<f64, 8, 8> {
        let a: Matrix4<f64> = 0.5
            * Matrix4::new(
                -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0,
            );
        let id = Matrix2::identity();
        let res = a.kronecker(&id);

        dbg!("amplifier:");
        println!("{}", res);
        res
    }

    fn grover() {
        let mut inputs = [0.0; 8].to_vec();
        inputs[1] = 1.0; // 001 -> 0100000000

        let mut state = EntangledParticleN::new(SVector::<f64, 8>::from_vec(inputs));

        let hadamard2x2 = SingleInputGate::Hadamard.get_matrix::<2>(0);
        let hadamard8x8: Matrix<f64, Const<8>, Const<8>, ArrayStorage<f64, 8, 8>> =
            hadamard2x2.kronecker(&hadamard2x2).kronecker(&hadamard2x2);

        // hadamard to input
        // puts inputs into superposition

        state.change_state_by_matrix(hadamard8x8);
        dbg!(&state);
        state.change_state_by_matrix(oracle());
        dbg!(&state);
        state.change_state_by_matrix(amplifier());
        dbg!(&state);

        state.measure(0);
        state.measure(1);
        state.measure(2);
        state.measure(3);
        state.measure(4);
        state.measure(5);
        state.measure(6);
        state.measure(7);
        dbg!(&state);
    }
    #[test]
    fn test_grover() {
        grover();
    }
}
