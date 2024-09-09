
/*
 * Benjamin Waters 
 * Nuclear Physics
 * Liquid Drop Model Binding Energy Calculato
 *
 * 9/5/2024
 *
*/

#![allow(non_upper_case_globals)]   // for scientific units
#![allow(non_snake_case)]           // for element names such as U for Uranium

use elements::Isotope;

pub mod elements;


fn hw1_isotopes() {
    
    let U_236 = Isotope::from_nucleons(236, 92);

    let Pd_117 = Isotope::from_nucleons(117, 46);

    let Xe_140 = Isotope::from_nucleons(140, 54);

    let Sr_94 = Isotope::from_nucleons(94, 38);
    

    U_236.report();
    Pd_117.report();
    Xe_140.report();
    Sr_94.report();   
}



fn main(){
    
    hw1_isotopes();


}


