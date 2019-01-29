use std::ops::{Add, Sub, Mul, Div, Neg};
use std::fmt;

use crate::finite_fields::FieldValues;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct KPoint<K> {
    pub x: K,
    pub y: K,
}

impl<K> fmt::Display for KPoint<K>
    where K : fmt::Display{
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "({}, {})", self.x, self.y)
        }
}


#[derive(Debug, PartialEq, Clone, Copy)]
pub enum ProjKPoint<K> {
    FinPoint(KPoint<K>),
    InfPoint,
}

impl<K> fmt::Display for ProjKPoint<K>
    where K : fmt::Display{
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            use ProjKPoint::*;
            match self {
                FinPoint(p) => write!(f, "FP: {}", p),
                InfPoint => write!(f, "Inf P"),
            }
        }
}

#[derive(Debug, PartialEq)]
pub struct EllipticCurve<K> {
    a: K,
    b: K
}

impl<K> fmt::Display for EllipticCurve<K>
    where K : fmt::Display{
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "y² = x^3 + {}*x + {}", self.a, self.b)
        }
}

impl<K> EllipticCurve<K> 
    where K : Add<Output=K> + Sub<Output=K> + Mul<Output=K> + Div<Output=K> + Neg<Output=K>
     + FieldValues<K> + PartialEq + Copy + fmt::Display{
        pub fn new(a: K, b: K) -> EllipticCurve<K>{
            EllipticCurve::<K>{
                a, b
            }
        }

        pub fn is_on_curve(&self, point : &ProjKPoint<K>) -> bool{
            use ProjKPoint::*;
            match point{
                InfPoint => true,
                FinPoint(p) =>{
                    let x = p.x;
                    let y = p.y;
                    y*y == x*x*x + self.a*x + self.b
                }
            }
        }

        pub fn neg_point(&self, point : ProjKPoint<K>) -> ProjKPoint<K>{
            assert!(self.is_on_curve(&point));
            use ProjKPoint::*;
            match point {
                InfPoint => InfPoint,
                FinPoint(p) => FinPoint(KPoint::<K>{
                    x: p.x, y: -p.y,
                }),
            }
        }

        pub fn add_points(&self, point1 : ProjKPoint<K>, point2 : ProjKPoint<K>) -> ProjKPoint<K>{
            assert!(self.is_on_curve(&point1));
            assert!(self.is_on_curve(&point1));

            use ProjKPoint::*;
            let a = self.a;
            let b = self.b;
            if point1 == self.neg_point(point2){
                return InfPoint;
            }
            match (&point1, &point2){
                (InfPoint, _) => point2,
                (_, InfPoint) => point1,
                (FinPoint(p1), FinPoint(p2)) => {
                    let (lambda, nu) = 
                        if p1.x != p2.x {
                            ((p2.y-p1.y)/(p2.x-p1.x), (p1.y*p2.x - p2.y*p1.x)/(p2.x-p1.x))
                        }else{
                            ((K::from_int(3)*p1.x*p1.x + a)/(K::from_int(2)*p1.y), (-p1.x*p1.x*p1.x + a*p1.x + K::from_int(2)*b)/(K::from_int(2)*p1.y))
                        };
                    let x3 = lambda*lambda - p1.x - p2.x;
                    let y3 = -lambda*x3 - nu;
                    FinPoint(KPoint::<K>{
                        x: x3, y: y3
                    })
                }
            }
        }

        pub fn scalar_mult(&self, n : i32, point : ProjKPoint<K>) -> ProjKPoint<K>{
            use ProjKPoint::*;
            let mut m = n;
            let mut p = point;
            let mut p1 = self.add_points(p, p);

            assert!(self.is_on_curve(&point));
            if n == 0{
                return InfPoint;
            }
            if n < 0{
                m = -m;
                p = self.neg_point(point);
            }

            let mut logm = 0;
            let mut m1 = m;
            while m1 != 0{
                m1 >>= 1;
                logm += 1;
            }

            println!("{}, {:b}", logm, m);
            while logm != 0{
                let bit = (m&(1<<(logm-1)))>>(logm-1); // the current bit
                println!("{}", bit);
                println!("{}", p);
                logm -= 1;
                if bit == 0{
                    p1 = self.add_points(p, p1);
                    p = self.add_points(p, p);
                }else{
                    p = self.add_points(p, p1);
                    p1 = self.add_points(p1, p1);
                }
            }
            p
        }
}

#[cfg(test)]
mod test;