#![feature(type_alias_enum_variants)]

extern crate rand;

pub mod finite_fields;
pub mod elliptic_curves;

use finite_fields::*;

declare_finite_field!(F2, 2, m2);

declare_elliptic_curve!(E, 1, 18);

type Point = PointEllipticCurve<E>;

fn main() {
    let mut p = Point::new(GL::new(0), GL::new(1));
    println!("Finding a point...\n");
    while !(E::is_on_curve(&p)){
        p = Point::new(GL::new(rand::random::<finite_fields::Integer>()), 
                       GL::new(rand::random::<finite_fields::Integer>()));
    }
    println!("Found ! computing order...\n");
    let mut p1 = p;

    let mut i = 1;

    while p1 != Point::InfPoint && E::is_on_curve(&p1){
        println!("{}", p1);
        p1 = p1 + p;
        i+=1;
    }

    if !E::is_on_curve(&p1){
        println!("Error");
    }
    println!("\nOrder of {} : {}", p, i);
}
