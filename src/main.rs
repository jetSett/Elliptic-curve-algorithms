#![feature(type_alias_enum_variants)]

extern crate rand;

pub mod finite_fields;
mod elliptic_curves;

use elliptic_curves::*;

#[derive(Debug)]
struct E{
}

type Point = PointEllipticCurve<E>;

impl EllipticCurve<E> for E {
    fn get_a_weierstrass() -> GL{
        GL::new(1)
    }
    fn get_b_weierstrass() -> GL{
        GL::new(1841)
    }
}


fn main() {
    let mut p = Point::new(GL::new(0), GL::new(1));
    println!("Finding a point...\n");
    while !(E::is_on_curve(&p)){
        p = Point::new(GL::new(rand::random::<i32>()), GL::new(rand::random::<i32>()));
    }
    println!("Found ! computing order...\n");
    let mut p1 = p;
    let mut p1_old = p;

    let mut i = 1;

    while p1 != Point::InfPoint && E::is_on_curve(&p1){
        p1_old = p1;
        p1 = p1 + p;
        println!("{}", E::is_on_curve(&p1));
        println!("{}", p1);
        i+=1;
    }

    println!("{}", p1_old);

    if! E::is_on_curve(&p1){
        println!("Error");
    }
    println!("{}", i);
}
