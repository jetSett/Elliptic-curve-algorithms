use std::marker::PhantomData;

use std::clone::Clone;

use std::ops::{Add, Sub, Neg};


use super::finite_fields::{Fp, IntegerAsType, Integer};

#[derive(Debug)]
pub struct TypeInt{
}

impl IntegerAsType for TypeInt{
    fn value() -> Integer{
        6569
    }
}

pub type GL = Fp<TypeInt>;

#[derive(Debug)]
pub struct GLPointEC<E : EllipticCurve<E>>{
    x : GL,
    y : GL,

    _phantom : PhantomData<E>
}

impl<E> GLPointEC<E>
    where E : EllipticCurve<E>{
        pub fn new(x : GL, y : GL) -> GLPointEC<E>{
            GLPointEC::<E>{
                x,
                y,

                _phantom: PhantomData
            }
        }
}

impl<E> PartialEq for GLPointEC<E>
    where E : EllipticCurve<E>{
        fn eq(&self, other: &GLPointEC<E>) -> bool{
            self.x == other.x && self.y == other.y
        }
}


impl<E> Clone for GLPointEC<E>
    where E : EllipticCurve<E>{
        fn clone(&self) -> GLPointEC<E>{
            GLPointEC::new(self.x, self.y)
        }
}

impl<E> Copy for GLPointEC<E>
    where E : EllipticCurve<E>{
}

#[derive(Debug)]
pub enum PointEllipticCurve<E : EllipticCurve<E>>{
    FinPoint(GLPointEC<E>),
    InfPoint,
}

use PointEllipticCurve::{FinPoint, InfPoint};

impl<E> PointEllipticCurve<E>
    where E : EllipticCurve<E>{
        pub fn new(x : GL, y : GL) -> PointEllipticCurve<E>{
            PointEllipticCurve::<E>::FinPoint(GLPointEC{
                x,
                y,

                _phantom: PhantomData
            })
        }

        pub fn infinite() -> PointEllipticCurve<E>{
            PointEllipticCurve::<E>::InfPoint
        }
}

impl<E> Clone for PointEllipticCurve<E>
    where E : EllipticCurve<E>{
        fn clone(&self) -> PointEllipticCurve<E>{
            match self{
                InfPoint => InfPoint,
                FinPoint(p) => FinPoint(*p)
            }
        }
}

impl<E> Copy for PointEllipticCurve<E>
    where E : EllipticCurve<E>{

}

pub trait EllipticCurve<E : EllipticCurve<E>>{

    fn get_a_weierstrass() -> GL;
    fn get_b_weierstrass() -> GL;

    fn is_on_curve(point : &PointEllipticCurve<E>) -> bool{
        match point {
            PointEllipticCurve::FinPoint(p) => {
                let x = p.x;
                let y = p.y;

                let a = E::get_a_weierstrass();
                let b = E::get_b_weierstrass();

                y*y == x*x*x + a*x + b
            },
            PointEllipticCurve::InfPoint => true,
        }
    }

    fn add_points(point1 : PointEllipticCurve<E>, point2 : PointEllipticCurve<E>) -> PointEllipticCurve<E>{
        if point2 == E::neg_point(&point1){
            println!("Coucou\n");
            return InfPoint;
        }

        let a = E::get_a_weierstrass();
        let b = E::get_b_weierstrass();

        match (&point1, &point2){
            (InfPoint, _) => point2,
            (_, InfPoint) => point1,
            (FinPoint(p1), FinPoint(p2)) => {
                let (lambda, nu) = 
                    if p1.x != p2.x {
                        ((p2.y-p1.y)/(p2.x-p1.x), (p1.y*p2.x - p2.y*p1.x)/(p2.x-p1.x))
                    }else{
                        ((GL::new(3)*p1.x*p1.x)/(GL::new(2)*p1.y), (-p1.x*p1.x*p1.x + a*p1.x + GL::new(2)*b)/(GL::new(2)*p1.y))
                    };
                let x3 = lambda*lambda - p1.x - p2.x;
                let y3 = -lambda*x3 - nu;
                FinPoint(GLPointEC::new(x3,y3))
            }
        }
    }

    fn neg_point(point : &PointEllipticCurve<E>) -> PointEllipticCurve<E>{
        match point {
            InfPoint => InfPoint,
            FinPoint(p) => FinPoint(GLPointEC::new(p.x, -p.y)),
        }
    }
}

impl<E> PartialEq for PointEllipticCurve<E>
    where E : EllipticCurve<E>{
        fn eq(&self, other: &PointEllipticCurve<E>) -> bool{
            match (self, other){
                (FinPoint(p1), FinPoint(p2)) => p1 == p2,
                (InfPoint, InfPoint) => true,
                _ => false
            }
        }
}

impl<E> Add for PointEllipticCurve<E>
    where E : EllipticCurve<E>{
        type Output = PointEllipticCurve<E>;

        fn add(self, other: PointEllipticCurve<E>) -> PointEllipticCurve<E>{
            E::add_points(self, other)
        }
}

impl<E> Neg for PointEllipticCurve<E>
    where E : EllipticCurve<E>{
        type Output = PointEllipticCurve<E>;

        fn neg(self) -> PointEllipticCurve<E>{
            E::neg_point(&self)
        }
}

impl<E> Sub for PointEllipticCurve<E>
    where E : EllipticCurve<E>{
        type Output = PointEllipticCurve<E>;

        fn sub(self, other: PointEllipticCurve<E>) -> PointEllipticCurve<E>{
            E::add_points(self, -other)
        }
}

