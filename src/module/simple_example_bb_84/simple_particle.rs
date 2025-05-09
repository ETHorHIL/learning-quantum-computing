use super::utils::round_to_n_decimal_places;
use nalgebra::{DVector, Matrix2, Matrix4, SVector, Vector2, Vector4};
use rand::prelude::*;
use std::{f64::consts::PI, vec};

#[derive(Default)]
pub struct System {
    pub apparatus: Apparatus,
    pub particle: Particle,
}

impl System {
    pub fn new(apparatus: Apparatus, particle: Particle) -> Self {
        Self {
            apparatus,
            particle,
        }
    }
    pub fn measure(&mut self) -> bool {
        // if the particle has not been measured before and we dont know its state we measure it now
        // orientation is true for north or false for south
        // state vector is the first (north) or second (south) basis vector of the appratus
        if self.particle.state_vector.is_none() {
            // 50% chance of 1.0, 50% chance of 0.0
            return self.set_particle_orientation_and_state(0.5);
        }
        // if the particles state vector is not one of the basis vectors we return a new state vector based on its probability
        self.set_particle_orientation_and_state(self.probability_north())
    }

    fn set_particle_orientation_and_state(&mut self, probability_north: f64) -> bool {
        let orientation = rand::rng().random_bool(probability_north);
        self.particle.orientation = Some(orientation);

        if orientation {
            self.particle.state_vector = Some(self.apparatus.basis.column(0).into());
        } else {
            self.particle.state_vector = Some(self.apparatus.basis.column(1).into());
        }
        orientation
    }

    pub fn probability_north(&self) -> f64 {
        // first get the vector basis * state
        let probability_north =
            self.apparatus.basis.transpose() * self.particle.state_vector.unwrap();
        //take square of first element to get probability for north
        probability_north[0].powi(2)
    }
}

#[derive(Default, Clone, Debug)]
pub struct Particle {
    state_vector: Option<Vector2<f64>>,
    orientation: Option<bool>,
}

pub struct Apparatus {
    basis: Matrix2<f64>,
    angle: f64,
}

impl Default for Apparatus {
    fn default() -> Apparatus {
        Apparatus {
            basis: Matrix2::new(1.0, 0.0, 0.0, 1.0),
            angle: 0.0,
        }
    }
}

impl Apparatus {
    /// rotate by teta 180degree should result in opposite direction
    pub fn set_angle(&mut self, teta: f64) {
        // convert teta angle to radians
        // PI / 2.0 is 90 degres in radian
        self.angle = teta;
        let mut radians = 0.0;
        if teta != 0.0 {
            radians = PI / (180.0 / (teta / 2.0));
        }

        self.basis = Matrix2::new(
            radians.cos(),
            radians.sin(),
            -1.0 * radians.sin(),
            radians.cos(),
        );
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_apparatus_rotation() {
        // apparatus
        let mut app = Apparatus::default();
        println!("basis {}", app.basis);
        app.set_angle(180.0);

        println!("basis{}", app.basis); // should be 0 1 -1 0
        app.set_angle(0.0);
        println!("basis{}", app.basis); // should be 1 0 1 0
    }

    #[test]
    fn test_measurement() {
        let apparatus = Apparatus::default();
        let particle = Particle::default();
        let mut system = System::new(apparatus, particle.clone());

        let obs_0 = system.measure();

        println!("first measurement {obs_0}");
        for _i in 0..10 {
            let obs_n = system.measure();
            assert_eq!(obs_0, obs_n);
        }
        println!("repeating previous measurement 10 times gave same answer");
        println!("rotate apparatus by 180 degree expect opposite outcome");

        system.apparatus.set_angle(180.0);
        let obs_1 = system.measure();
        assert_eq!(obs_0, !obs_1);
        println!("first measurement is indeed opposite outcome {obs_0}");
        for _i in 0..10 {
            let obs_n = system.measure();
            assert_eq!(obs_1, obs_n);
        }
        println!("repeating measurement ten times give same answer");

        // 1.0 at 0.0 gives true as does 0.70,-0.7

        let mut res = false;
        system.apparatus.set_angle(0.0);
        res = system.measure();
        let mut counter = 0;
        while !res {
            system.apparatus.set_angle(90.0);
            system.measure();
            system.apparatus.set_angle(0.0);
            res = system.measure();
            counter += 1;
            assert_ne!(counter, 10);
        }
        for _i in 0..9 {
            assert!(system.measure());
        }

        res = false;
        system.apparatus.set_angle(90.0);
        res = system.measure();

        let mut counter = 0;
        while !res {
            system.apparatus.set_angle(0.0);
            system.measure();
            system.apparatus.set_angle(90.0);
            res = system.measure();
            counter += 1;
            assert_ne!(counter, 10);
        }
        for i in 0..9 {
            if !system.measure() {
                panic!("");
            }
        }
        dbg!("1.0 gives at degree true as does 0.7,-0.7 at 90");

        // reset apparatus and particle state
        system.apparatus.set_angle(0.0);
        system.measure();

        let rotation_angle = 45.0;

        println!("rotate apparatus by {rotation_angle} degree after each measurement, expect equal probability for same outcome");
        let mut trues = 0;
        let mut falses = 0;

        for _i in 0..10 {
            let new_angle = (system.apparatus.angle + rotation_angle).rem_euclid(360.0);

            system.apparatus.set_angle(new_angle);
            let obs_2 = system.measure();
            if obs_2 {
                trues += 1;
            } else {
                falses += 1;
            }
        }

        println!("after 100 tests we have {trues} times N and {falses} times S");
    }
}
