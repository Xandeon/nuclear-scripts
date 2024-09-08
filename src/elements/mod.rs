

/*
 * Benjamin Waters
 * Nuclear Physics
 *
 * Elements Module
 * 9/8/2024
 *
 */

#![allow(non_upper_case_globals)]   // for scientific units
#![allow(non_snake_case)]           // for element names such as U for Uranium

// ---------------------------- units ----------------------------------

const eV: f64 = 1.0;
const keV: f64 = 1000.0;
const MeV: f64 = 1000.0 * 1000.0;
const GeV: f64 = 1000.0 * 1000.0 * 1000.0;

// ---------------------------- nucleons constants ---------------------

const m_electron: f64 = 0.511; // MeV/c^2
const m_proton: f64 = 938.280; // MeV/c^2q
const m_neutron: f64 = 939.573; // MeV/c^2

const amu: f64 = 931.5; // 1 amu = 931.5 MeV/c^2 = 1.661 * 10e-27 kg

// ---------------------------- liquid drop Model constants ------------

const A_V: f64 = 15.5; // MeV (volume coefficient)
const A_S: f64 = 16.8; // MeV (surface area coeffiecient)
const A_C: f64 = 0.72; // eV (coulomb coefficient)
const A_A: f64 = 23.0; // MeV (asymmetry term coefficient)
const A_P: f64 = 34.0; // MeV (pairing term coefficient)    


pub struct Isotope{
    element: String,
    A: i32,
    Z: i32,
    N: i32,
}

impl Isotope{
    pub fn create_isotope(element: String, A: i32, Z: i32, N: i32) -> Self{
        Self{
            element, A, Z, N,
        }
    }
    
    pub fn get_upper_isobar(&self) -> Self{
        Isotope::create_isotope(String::from("UnIm"), self.A, self.Z + 1, self.N - 1) 
    }

    pub fn mass(&self) -> f64{
        // returns the mass in Mev/c^2 units, divide by 'amu' to get amu units
        (self.Z as f64) * m_proton +
        (self.N as f64) * m_neutron -
        self.binding_energy()
    }
    
    pub fn binding_energy(&self) -> f64{
        /// Liquid Drop Model
        /// Implements the Bethe-Weinsaecker semi-empirical mass formula for binding energy
        ///
        /// The calculation is performed in units of MeV and can be converted afterward to 
        /// other units to not clutter the calculation. Only the coloumb term is scaled.
        
        let Isotope {A, Z, N, ..} = self; 
        let A_f64: f64 = *A as f64;
        let Z_f64: f64 = *Z as f64;
        let N_f64: f64 = *N as f64;

        let volume = A_V * A_f64;
        let surface = -1.0 * A_S * A_f64.powf(2.0 / 3.0);
        let coulomb = -1.0 * A_C * Z_f64.powf(2.0) * A_f64.powf(-1.0 / 3.0); 
        let asymetry = -1.0 * A_A * (N_f64 - Z_f64).powf(2.0)/A_f64;
        let delta = A_P * A_f64.powf(-3.0 / 4.0);

        let pairing ={
            if Z % 2 == 0 && N % 2 == 0{            // Even Z / Even N
                delta
            } else if Z % 2 != 0 && N % 2 == 0{     // Odd Z / Even N
                0.0
            } else if Z % 2 != 0 && N % 2 != 0{     // Odd Z / Odd N
                -1.0 * delta
            } else {                                // Not a common case
                0.0    
            } 
        };
    
        volume + surface + coulomb + asymetry + pairing // return the binding energy
    }

    pub fn binding_energy_per_nucleon(&self) -> f64{
        // Returns the binding energy per nucleon of the Isotope
    
        self.binding_energy() / self.A as f64     
    }

    pub fn report(&self) {
        println!(
                "------------- Isotope Report ------------ \n\
                 {:<20} {:>15} \n\
                 {:<20} {:>15} nucleons \n\
                 {:<20} {:>15} protons \n\
                 {:<20} {:>15} neutrons \n\
                 {:<20} {:>15.5} amu \n\
                 {:<20} {:>15.5} MeV/c^2 \n\
                 
                 {:<20} {:>15.5} MeV \n\
                 {:<20} {:>15.5} Mev \n\
                 -----------------------------------------", 
                 "Element:", self.element,
                 "A", self.A,
                 "Z", self.Z,
                 "N", self.N,
                 "Mass (amu)", self.mass() / amu,
                 "", self.mass(),
                 "Binding Energy:", self.binding_energy(),
                 "BE/A:", self.binding_energy_per_nucleon(),
                ); 
    }

}



