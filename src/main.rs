#![feature(type_alias_enum_variants)]

#[macro_use] pub mod finite_fields;
pub mod elliptic_curves;

use finite_fields::*;
use elliptic_curves::*;

declare_finite_field!(K, 10169, m10169);

fn main() {

    let ell = EllipticCurve::<K>::new(K::new(1), K::new(5841));
    let mut p = ProjKPoint::FinPoint(KPoint::<K>{
        x: K::new(0), y: K::new(0)
    });

    while !(ell.is_on_curve(&p)) {
        p = ProjKPoint::FinPoint(KPoint{
            x: K::new(rand::random::<finite_fields::Integer>()),
            y: K::new(rand::random::<finite_fields::Integer>()),
        });
    }
    
    println!("Found ! computing order...\n");
    let mut p1 = p;

    let mut i = 1;

    while p1 != ProjKPoint::InfPoint && ell.is_on_curve(&p1){
        println!("{}", p1);
        p1 = ell.add_points(p1, p);
        i+=1;
    }

    if !ell.is_on_curve(&p1){
        println!("Error");
    }
    println!("\nOrder of {} : {}", p, i);
}
