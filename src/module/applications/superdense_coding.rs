use crate::module::entangled_particle_n::EntangledParticleN;
use crate::module::gates::basic::*;
use nalgebra::{ArrayStorage, Const, Matrix, Matrix2, Matrix4, SVector, Vector2, Vector4};

#[cfg(test)]
mod tests {

    use crate::module::{
        gates::{
            self,
            basic::{self, SingleInputGate},
        },
        utils::round_to_n_decimal_places,
    };

    use super::*;

    #[test]
    fn test_superdense_coding() {
        // Alice wants to send two classical bits of information using only one qbit
        // variants of two classical bits
        // two entangled particles in bell stared are prepared
        // one is given to alice and one to bob
        let prtcl_template = EntangledParticleN::new(SVector::<f64, 4>::new(
            0.5f64.sqrt(),
            0.0,
            0.0,
            0.5f64.sqrt(),
        ));

        let mut prtcl_00 = prtcl_template.clone();

        // gates
        let hadamard = SingleInputGate::Hadamard.get_matrix::<4>(0);
        let cnot = TwoInputGate::CNot.get_matrix::<4>(0, 1);
        let x_gate = SingleInputGate::X.get_matrix::<4>(0);
        let z_gate = SingleInputGate::Z.get_matrix::<4>(0);
        let y_gate = SingleInputGate::Y.get_matrix::<4>(0);

        // alice applies circuit to her qbit depending on which bit pair she wants to send

        // in the 00 case alice does nothing to her qbits
        // in the 01 case alice applies the x gate
        // for 10 she applies the z gate
        // for 11 she applies the y gate

        //00 do nothing
        // bob applies reverse bell circuit
        prtcl_00.change_state_by_matrix(cnot);
        prtcl_00.change_state_by_matrix(hadamard);
        // bob measures alices qubit
        let res_0 = prtcl_00.measure(0);
        // bob measures his own qubit
        let res_1 = prtcl_00.measure(1);
        // the result should be 0,0
        assert_eq!([res_0, res_1], [true, true]);

        // 01
        // prepare new particle
        let mut prtcl_01 = prtcl_template.clone();

        // alice applies the X gate
        prtcl_01.change_state_by_matrix(x_gate);

        // bob applies reverse bell circuit
        prtcl_01.change_state_by_matrix(cnot);
        prtcl_01.change_state_by_matrix(hadamard);
        // bob measures alices qubit
        let res_0 = prtcl_01.measure(0);
        // bob measures his own qubit
        let res_1 = prtcl_01.measure(1);
        // the result should be 1,0
        assert_eq!([res_0, res_1], [true, false]);

        // 10
        // prepare new particle
        let mut prtcl_10 = prtcl_template.clone();
        // alice applies the X gate
        prtcl_10.change_state_by_matrix(z_gate);

        // bob applies reverse bell circuit
        prtcl_10.change_state_by_matrix(cnot);
        prtcl_10.change_state_by_matrix(hadamard);
        // bob measures alices qubit
        let res_0 = prtcl_10.measure(0);
        // bob measures his own qubit
        let res_1 = prtcl_10.measure(1);
        // the result should be 0,1
        assert_eq!([res_0, res_1], [false, true]);

        // 11
        // prepare new particle
        let mut prtcl_11 = prtcl_template.clone();
        // alice applies the X gate
        prtcl_11.change_state_by_matrix(y_gate);

        // bob applies reverse bell circuit
        prtcl_11.change_state_by_matrix(cnot);
        prtcl_11.change_state_by_matrix(hadamard);

        // bob measures alices qubit
        let res_0 = prtcl_11.measure(0);
        // bob measures his own qubit
        let res_1 = prtcl_11.measure(1);
        // the result should be 1,1
        assert_eq!([res_0, res_1], [false, false]);
        dbg!("worked");
    }
}
