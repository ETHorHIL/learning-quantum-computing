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

mod tests {

    use super::*;

    #[test]
    fn test_teleport_pure_0() {
        // Alice wants to send two classical bits of information using only one qbit
        // variants of two classical bits
        // two entangled particles in bell stared are prepared
        // one is given to alice and one to bob

        // alice and bob share a particle in bell state
        // alice also has a particle in purestate 0 (1,0)
        let initial_state: SVector<f64, 8> = SVector::from_vec(vec![
            0.5f64.sqrt(),
            0.0,
            0.0,
            0.5f64.sqrt(),
            0.0,
            0.0,
            0.0,
            0.0,
        ]);
        let mut prtcl = EntangledParticleN::new(initial_state);

        // alice applies the reverse bell circuit and measures
        prtcl.change_state_by_matrix(TwoInputGate::CNot.get_matrix::<8>(0, 1));

        prtcl.change_state_by_matrix(SingleInputGate::Hadamard.get_matrix::<8>(0));

        let alice_measurement = [prtcl.measure(0), prtcl.measure(1)];

        // alice sends her measurement result to bob

        // in the 00 case bob knows his qbits have been teleported correctly and does nothing her qbits
        // in the 01 case bob applies the x gate
        // for 10 bob applies the z gate
        // for 11 bob applies the y gate
        match alice_measurement {
            [true, true] => {}
            [true, false] => prtcl.change_state_by_matrix(SingleInputGate::X.get_matrix::<8>(2)),
            [false, true] => prtcl.change_state_by_matrix(SingleInputGate::Z.get_matrix::<8>(2)),
            [false, false] => prtcl.change_state_by_matrix(SingleInputGate::Y.get_matrix::<8>(2)),
        }

        let params = prtcl.get_params();
        assert_eq!(round_to_n_decimal_places(params.sum(), 5), 1.0);
        for (index, amplitude) in params.iter().enumerate() {
            if index % 2 == 0 {
                assert!(
                    round_to_n_decimal_places(*amplitude, 5) == 0.0
                        || round_to_n_decimal_places(*amplitude, 5) == 1.0
                );
            } else {
                assert!(round_to_n_decimal_places(*amplitude, 5) != 1.0);
            }
        }
    }

    #[test]
    fn test_teleport_superposition() {
        // Alice wants to send two classical bits of information using only one qbit
        // variants of two classical bits
        // two entangled particles in bell stared are prepared
        // one is given to alice and one to bob

        // alice and bob share a particle in bell state
        // alice also has a particle in superposition 1/sqrt(2) * (1,0)+1/sqrt(2) * (0,1)
        let initial_state: SVector<f64, 8> = SVector::from_vec(vec![
            0.5f64.sqrt(),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.5f64.sqrt(),
        ]);
        let mut prtcl = EntangledParticleN::new(initial_state);

        // alice applies the reverse bell circuit and measures
        prtcl.change_state_by_matrix(TwoInputGate::CNot.get_matrix::<8>(0, 1));
        prtcl.change_state_by_matrix(SingleInputGate::Hadamard.get_matrix::<8>(0));

        let alice_measurement = [prtcl.measure(0), prtcl.measure(1)];
        dbg!(&prtcl);
        // alice sends her measurement result to bob

        // in the 00 case bob knows his qbits have been teleported correctly and does nothing her qbits
        // in the 01 case bob applies the x gate
        // for 10 bob applies the z gate
        // for 11 bob applies the y gate
        match alice_measurement {
            [true, true] => {}
            [true, false] => prtcl.change_state_by_matrix(SingleInputGate::X.get_matrix::<8>(2)),
            [false, true] => prtcl.change_state_by_matrix(SingleInputGate::Z.get_matrix::<8>(2)),
            [false, false] => prtcl.change_state_by_matrix(SingleInputGate::Y.get_matrix::<8>(2)),
        }

        dbg!(&prtcl);

        let params = prtcl.get_params();
        assert_eq!(
            round_to_n_decimal_places((0.5f64.sqrt() * params).sum(), 5),
            1.0
        );
        for amplitude in params {
            assert!(
                round_to_n_decimal_places(*amplitude, 5) == 0.0
                    || round_to_n_decimal_places(*amplitude, 5)
                        == round_to_n_decimal_places(0.5f64.sqrt(), 5)
            );
        }
    }
}
