use std::marker::PhantomData;

use super::finite_fields::{Fp, IntegerAsType, Integer};

pub struct TypeInt{
}

impl IntegerAsType for TypeInt{
    fn value() -> Integer{
        29
    }
}

pub type GL = Fp<TypeInt>;

pub struct PointEllipticCurve<E : EllipticCurve<E>>{
    x : GL,
    y : GL,

    _phantom : PhantomData<E>
}

impl<E> PointEllipticCurve<E>
    where E : EllipticCurve<E>{
        pub fn new(x : GL, y : GL) -> PointEllipticCurve<E>{
            PointEllipticCurve::<E>{
                x,
                y,

                _phantom: PhantomData
            }
        }
}

pub trait EllipticCurve<E : EllipticCurve<E>>{

    fn get_a_weierstrass() -> GL;
    fn get_b_weierstrass() -> GL;

    fn is_on_curve(p : &PointEllipticCurve<E>) -> bool{
        let x = &p.x;
        let y = &p.y;
        y*y == &(x*x)*x + &E::get_a_weierstrass()*x + E::get_b_weierstrass()
    }

    fn add_points(p1 : &PointEllipticCurve<E>, p2 : &PointEllipticCurve<E>) -> PointEllipticCurve<E>{
        let a = E::get_a_weierstrass();
        let b = E::get_b_weierstrass();
    }

    fn neg_point(p1 : &PointEllipticCurve<E>){

    }

    fn 
}