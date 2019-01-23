#![feature(type_alias_enum_variants)]

pub mod finite_fields;
mod elliptic_curves;

use elliptic_curves::*;

#[derive(Debug)]
struct E11{
}

type Point = PointEllipticCurve<E11>;

impl EllipticCurve<E11> for E11 {
    fn get_a_weierstrass() -> GL{
        GL::new(1)
    }
    fn get_b_weierstrass() -> GL{
        GL::new(1)
    }
}


fn main() {
    let mut p3 = Point::new(GL::new(0), GL::new(1));

    println!("{}", E11::is_on_curve(&p3));

    let mut i = 0;

    assert!(GL::new(-1) == GL::new(4));

    while p3 != Point::InfPoint{
        p3 = p3 + p3;
        println!("{}", E11::is_on_curve(&p3));
        println!("{:?}", p3);
        i+=1;
    }
}
