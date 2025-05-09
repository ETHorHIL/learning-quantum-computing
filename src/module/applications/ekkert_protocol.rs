use super::utils::round_to_n_decimal_places;
use crate::module::entangled_particle_n::{self, EntangledParticleN};
use nalgebra::{DVector, Matrix2, Matrix4, SVector, Vector2, Vector4};
use rand::prelude::*;
use std::{f64::consts::PI, vec};

#[derive(Debug)]
pub struct EntangledParticleStream {
    particles: Vec<EntangledParticleN<4>>,
}

impl EntangledParticleStream {
    fn new(len: usize) -> Self {
        let mut particles = vec![];
        for i in 0..len {
            let prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(
                0.5f64.sqrt(),
                0.0,
                0.0,
                0.5f64.sqrt(),
            ));
            particles.push(prtcl);
        }
        Self { particles }
    }
    fn rotate_to_basis(basis: &Basis) -> f64 {
        match basis {
            Basis::Standard => 0.0,
            Basis::Degree90 => 90.0,
            Basis::Degree45 => 45.0,
        }
    }

    fn rotate_bases(&mut self, basis: &Vec<Basis>, alice: bool) {
        for i in 0..basis.len() {
            self.particles[i].swap_basis(basis[i].degrees(), 0);
        }
    }

    fn measure_all(&mut self, alice: bool) -> Vec<bool> {
        let mut res: Vec<bool> = vec![];

        for particle in &mut self.particles {
            res.push(particle.measure(0));
        }

        res
    }
}

#[derive(Debug)]
pub struct Party {
    name: String,
    bases: Option<Vec<Basis>>,
    measurements: Option<Vec<bool>>,
    key: Option<Vec<bool>>,
    non_agreemets: Option<Vec<bool>>,
}

impl Party {
    pub fn new(name: &str, no_particles: usize) -> Self {
        let mut party = Party {
            name: name.to_string(),
            bases: None,
            measurements: None,
            key: None,
            non_agreemets: None,
        };
        party.generate_random_bases_vec(no_particles);
        party
    }
    fn random_basis(&self) -> Basis {
        let x = rand::rng().random_range(0..3);
        match x {
            0 => Basis::Standard,
            1 => Basis::Degree90,
            _ => Basis::Degree45,
        }
    }

    fn generate_random_bases_vec(&mut self, len: usize) {
        let mut vec: Vec<Basis> = vec![];
        for _i in 0..len {
            vec.push(self.random_basis());
        }
        self.bases = Some(vec);
    }

    fn record_measurement(&mut self, obs: &Vec<bool>) {
        self.measurements = Some(obs.to_owned());
    }

    /// compare other basis to ones own. store measurements where bases agreed
    /// Store measurements where bases didnt agree in non_agreements
    fn compare_basis(&mut self, other: Vec<Basis>) {
        let mut key: Vec<bool> = vec![];
        let measurements = self.measurements.as_ref().expect("no measurements");
        let mut non_agreemets: Vec<bool> = vec![];
        let bases = self.bases.as_ref().expect("no bases");

        for (i, basis) in bases.iter().enumerate() {
            // if the bases agree, we store the measurement. This will be the key
            if basis == &other[i] {
                key.push(measurements[i]);
            }
            // if the bases dont agree we store the measuremens.
            else {
                non_agreemets.push(measurements[i]);
            }
        }
        self.key = Some(key);
        self.non_agreemets = Some(non_agreemets);
    }

    /// compare the measurements that were taken when bases didnt agree. We expect that
    /// around 1/4th of these measurements should agree if noone evesdropped
    fn compare_non_agreements(&self, other: Vec<bool>) -> f64 {
        assert!(other.len() == self.non_agreemets.as_ref().unwrap().len());
        let mut counter = 0;
        let non_agreements = self
            .non_agreemets
            .as_ref()
            .expect("no non agreements found");
        for (i, item) in non_agreements.iter().enumerate() {
            if item == &other[i] {
                counter += 1;
            }
        }
        counter as f64 / other.len() as f64
    }

    fn share_bases(&self) -> Vec<Basis> {
        self.bases.clone().expect("no basis found")
    }
    fn share_non_agreements(&self) -> Vec<bool> {
        self.non_agreemets
            .clone()
            .expect("no non agreements found ")
    }
}

#[derive(PartialEq, Eq, Clone, Debug)]
enum Basis {
    Standard,
    Degree90,
    Degree45,
}

impl Basis {
    pub fn degrees(&self) -> f64 {
        match self {
            Basis::Standard => 0.0,
            Basis::Degree90 => 180.0 * (1.0 / 3.0),
            Basis::Degree45 => 180.0 * (2.0 / 3.0),
        }
    }
}

fn generate_particle_stream(len: usize) -> Vec<EntangledParticleN<4>> {
    let mut vec: Vec<EntangledParticleN<4>> = vec![];
    for i in 0..len {
        let prtcl = EntangledParticleN::new(SVector::<f64, 4>::new(
            0.5f64.sqrt(),
            0.0,
            0.0,
            0.5f64.sqrt(),
        ));
        vec.push(prtcl);
    }
    vec
}

// too make particle stream a struct
mod tests {

    use super::*;

    #[test]
    fn test_ekkert_no_eve() {
        /*
        current problem:
        - keys dont match
        - Something about the entanglement aint working right
         */

        let no_particles = 1000;

        // parties receive random basis at generations
        let mut alice = Party::new("alice", no_particles);
        let mut bob = Party::new("bob", no_particles);
        let mut prtcls = EntangledParticleStream::new(no_particles);

        // alice rotates her basis according to her random bases
        prtcls.rotate_bases(alice.bases.as_ref().unwrap(), true);
        alice.record_measurement(&prtcls.measure_all(true));
        prtcls.rotate_bases(bob.bases.as_ref().unwrap(), false);
        bob.record_measurement(&prtcls.measure_all(false));

        alice.compare_basis(bob.share_bases());
        bob.compare_basis(alice.share_bases());

        let alice_result = alice.compare_non_agreements(bob.share_non_agreements());

        println!("comparing results. both should be around 1/4");
        println!(
            "alice_result: {} of {} nonagreements and {} agreements",
            alice_result,
            alice.non_agreemets.unwrap().len(),
            alice.key.unwrap().len()
        );
    }

    #[test]
    fn test_ekkert_eve() {
        /*
        current problem:
        - keys dont match
        - Something about the entanglement aint working right
         */

        let no_particles = 1000;

        // parties receive random basis at generations
        let mut alice = Party::new("alice", no_particles);
        let mut bob = Party::new("bob", no_particles);
        let mut eve = Party::new("eve", no_particles);
        let mut prtcls = EntangledParticleStream::new(no_particles);

        // alice rotates her basis according to her random bases
        prtcls.rotate_bases(alice.bases.as_ref().unwrap(), true);
        alice.record_measurement(&prtcls.measure_all(true));

        // eve measures
        prtcls.rotate_bases(eve.bases.as_ref().unwrap(), false);
        eve.record_measurement(&prtcls.measure_all(false));

        // bob measures
        prtcls.rotate_bases(bob.bases.as_ref().unwrap(), false);
        bob.record_measurement(&prtcls.measure_all(false));

        alice.compare_basis(bob.share_bases());
        bob.compare_basis(alice.share_bases());

        let alice_result = alice.compare_non_agreements(bob.share_non_agreements());

        println!("comparing results. both should be around 3/8");
        println!(
            "alice_result: {} of {} nonagreements and {} agreements",
            alice_result,
            alice.non_agreemets.unwrap().len(),
            alice.key.unwrap().len()
        );
    }
}
