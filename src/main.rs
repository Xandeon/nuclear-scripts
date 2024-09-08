
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
use std::env;
use std::fs;    // filesystem


pub mod elements;


fn hw1_isotopes() {
    
    let U_236 = Isotope::create_isotope(String::from("Uranium"), 236, 92, 144);

    let Pd_117 = Isotope::create_isotope(String::from("Palladium"), 117, 46, 117 - 46);

    let Xe_140 = Isotope::create_isotope(String::from("Xenon"), 140, 54, 140 - 54);

    let Sr_94 = Isotope::create_isotope(String::from("Strontium"), 94, 38, 94 - 38);
    

    U_236.report();
    Pd_117.report();
    Xe_140.report();
    Sr_94.report();   
}



fn main(){
    



}


