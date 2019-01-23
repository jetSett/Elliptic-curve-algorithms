pub mod finite_fields;
mod elliptic_curves;

use elliptic_curves::*;

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

    fn add_points(p1 : &Point, p2 : &Point) -> Point{
        Point::new(GL::new(0), GL::new(1))
    }
}


fn main() {
    let p1 = Point::new(GL::new(1), GL::new(0));
    let p2 = Point::new(GL::new(1), GL::new(1));
    let p3 = Point::new(GL::new(0), GL::new(1));

    println!("{}", E11::is_on_curve(&p1));
    println!("{}", E11::is_on_curve(&p2));
    println!("{}", E11::is_on_curve(&p3));
}
