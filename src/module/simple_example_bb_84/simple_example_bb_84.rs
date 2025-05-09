// BB84 protocol for secure communication using quatum mechanical properties
use super::simple_particle::{self as pa, Particle};
use rand::prelude::*;

#[derive(Default)]
pub struct Party {
    name: String,
    message: Option<Vec<bool>>,
    basis_vec: Option<Vec<bool>>,
    basis_agreements: Vec<bool>,
    msg_agreements: Vec<bool>,
}

impl Party {
    pub fn new(name: &str) -> Self {
        Party {
            name: name.to_string(),
            ..Default::default()
        }
    }
    pub fn compare_msg_bits(&self, other_bits: Vec<bool>) -> bool {
        let own_bits = self.msg_agreements[0..self.msg_agreements.len() / 2].to_owned();
        println!("Publicly comparing half of the remaining message (after removing the part where basis didnt agree)");
        println!("If noone listended in we should have identical output. We can then use the other half of the msg for a common key");
        dbg!(&own_bits);
        dbg!(&other_bits);
        self.msg_agreements[0..self.msg_agreements.len() / 2] == other_bits
    }

    pub fn compare_bases(&mut self, other_basis: Vec<bool>) {
        let basis_vec = self.basis_vec.as_ref().expect("no basis");
        for (i, item) in other_basis.iter().enumerate() {
            if &basis_vec[i] == item {
                self.basis_agreements.push(other_basis[i]);
                self.msg_agreements.push(self.message.as_ref().unwrap()[i]);
            }
        }
        println!(
            "{} compares received basis with own and finds {} agreements. Around {} should be expected",
            self.name,
            self.msg_agreements.len(),basis_vec.len() / 2
        );
    }

    pub fn share_n_bits(&self) -> Vec<bool> {
        self.msg_agreements[0..self.msg_agreements.len() / 2].to_vec()
    }
    pub fn share_basis(&self) -> Vec<bool> {
        self.basis_vec.clone().unwrap()
    }

    pub fn read_qbits(&mut self, qbits: Vec<pa::Particle>) {
        println!("{} reads qbits", self.name);
        let mut system = pa::System::default();
        self.create_basis_vec(qbits.len());
        let basis_vec = self.basis_vec.as_ref().unwrap();
        let mut message: Vec<bool> = vec![];
        for i in 0..qbits.len() {
            let angle = if basis_vec[i] { 0.0 } else { 90.0 };
            system.apparatus.set_angle(angle);
            system.particle = qbits[i].clone();
            let res = system.measure();
            message.push(res);
        }
        self.message = Some(message);
    }

    pub fn set_message_and_bases(&mut self, message: Vec<bool>) {
        let len = message.len();
        if len % 4 != 0 {
            panic!("message lenght needs to be divisible by 4")
        }
        self.message = Some(message);
        self.create_basis_vec(len);
    }
    pub fn generate_particle_stream(&mut self) -> Vec<Particle> {
        let message = self.message.as_ref().expect("no msg").clone();
        let basis_vec = self.basis_vec.as_ref().expect("no basis").clone();
        let mut particle_stream: Vec<Particle> = vec![];
        println!("{} is generating {} particles", self.name, message.len());
        for i in 0..message.len() {
            particle_stream.push(self.generate_particle(message[i], basis_vec[i]));
        }
        println!("new particle stream is generated");

        particle_stream
    }
    fn generate_particle(&self, bit: bool, basis: bool) -> Particle {
        // select basis and orient apparatus accordingly
        let angle = if basis { 0.0 } else { 90.0 };
        let mut system = pa::System::default();
        system.apparatus.set_angle(angle);
        let mut set_state = system.measure();
        let mut particle = system.particle.clone();

        // measure particle until the bit we want to represent is set
        while set_state != bit {
            system = pa::System::default();
            system.apparatus.set_angle(angle);
            set_state = system.measure();
            particle = system.particle.clone();
        }

        particle
    }
    fn create_basis_vec(&mut self, len: usize) {
        let mut basis_vec = vec![];

        for _bit in 0..len {
            basis_vec.push(rand::rng().random_bool(0.5));
        }
        self.basis_vec = Some(basis_vec);
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_bb84_with_eve() {
        // eve intercepts communication
        let mut alice = Party::new("Alice");
        let mut bob = Party::new("Bob");
        let mut eve = Party::new("Eve");

        // alice decides what message to send she creates qbits that encode the message and shares publicly
        // she randomly assigns bases to encode the qubits
        // message lenght has to be 4n
        alice.set_message_and_bases([true; 100].to_vec());

        let qbits = alice.generate_particle_stream();

        eve.read_qbits(qbits);
        let qbits = eve.generate_particle_stream();

        // bob reads the message using random bases
        bob.read_qbits(qbits);

        // bob and alice share their bases publicly and compare. They keep only the bits where they accidentally used the same basis
        // they should have around 2n bits in common
        bob.compare_bases(alice.share_basis());
        alice.compare_bases(bob.share_basis());

        // bob and alice share half of the 2n bits of the message
        // since alice intercepted, there should be only n/2 agreements so comparing n bits should reveal eve listening in
        println!("expecting to fail since alice listened in");
        assert!(!bob.compare_msg_bits(alice.share_n_bits()));
    }

    #[test]
    fn test_bb84_no_eve() {
        // noone intercepts communication

        let mut alice = Party::new("Alice");
        let mut bob = Party::new("Bob");
        dbg!("start");

        // alice decides what message to send she creates qbits that encode the message and shares publicly
        // she randomly assigns bases to encode the qubits
        // message lenght has to be 4n
        alice.set_message_and_bases(vec![
            true, true, true, true, true, true, true, true, true, true, true, true, true, true,
            true, true, true, true, true, true,
        ]);

        dbg!("message set");
        let qbits = alice.generate_particle_stream();

        // bob reads the message using random bases
        bob.read_qbits(qbits);
        dbg!("qbits read");

        // bob and alice share their bases publicly and compare. They keep only the bits where they accidentally used the same basis
        // they should have around 2n bits in common
        bob.compare_bases(alice.share_basis());
        alice.compare_bases(bob.share_basis());
        dbg!("bases compared");

        // bob and alice share half of the 2n bits of the message
        // if noone intercepted there should be agreement
        assert!(bob.compare_msg_bits(alice.share_n_bits()));
    }
}
