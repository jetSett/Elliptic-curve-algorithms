use std::fmt;

use crate::field::{FieldValues, Field};

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum ProjKPoint<K> {
    FinPoint(K, K),
    InfPoint,
}

impl<K> fmt::Display for ProjKPoint<K>
    where K : fmt::Display{
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            use ProjKPoint::*;
            match self {
                FinPoint(x, y) => write!(f, "FP: ({}, {})", x, y),
                InfPoint => write!(f, "Inf P"),
            }
        }
}


// y^2 + a_1 xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6
#[derive(Debug, PartialEq)]
pub struct EllipticCurve<K> {
    pub a_1: K,
    pub a_3: K,

    pub a_2: K,
    pub a_4: K,
    pub a_6: K,
}

impl<K> EllipticCurve<K> 
    where K : Field
     + FieldValues<K>{
        pub fn new_reduced_weierstrass(a: K, b: K) -> EllipticCurve<K>{
            EllipticCurve::<K>{
                a_1: K::from_int(0),
                a_3: K::from_int(0),

                a_2: K::from_int(0),
                a_4: a,
                a_6: b,

            }
        }

        pub fn is_reduced_weierstrass(&self) -> bool{
            self.a_1 == K::from_int(0) && self.a_2 == K::from_int(0) && self.a_3 == K::from_int(0)
        }

        pub fn to_reduced_weierstrass(self) -> EllipticCurve<K>{
            EllipticCurve::new_reduced_weierstrass(-K::from_int(27)*self.c4(), -K::from_int(54)*self.c6())
        }

        fn b2(&self) -> K{
            self.a_1*self.a_1 + K::from_int(4)*self.a_2
            // self.a_1*self.a_1 + K::from_int(4)*self.a_4 // Version from Silvermann
        }

        fn b4(&self) -> K{
            K::from_int(2)*self.a_4 + self.a_1*self.a_3
        }

        fn b6(&self) -> K{
            self.a_3*self.a_3 + K::from_int(4)*self.a_6
        }

        fn c4(&self) -> K{
            let b2 = self.b2();
            let b4 = self.b4();

            b2*b2 - K::from_int(24)*b4
        }

        fn c6(&self) -> K{
            let b2 = self.b2();
            let b4 = self.b4();
            let b6 = self.b6();

            -b2*b2*b2 + K::from_int(36)*b2*b4 - K::from_int(216)*b6
        }

        pub fn j_invariant(&self) -> K{
            let c4 = self.c4();
            c4*c4*c4/self.discriminant()
        }

        pub fn discriminant(&self) -> K{
            let b2 = self.b2();
            let b4 = self.b4();
            let b6 = self.b6();

            K::from_int(-4)*(b2*b2*b2*b6 - b2*b2*b4*b4 - K::from_int(36)*b2*b4*b6 + K::from_int(32)*b4*b4*b4 + K::from_int(108)*b6*b6)
            //-b2*b2*b8 - K::from_int(8)*b4*b4*b4 - K::from_int(27)*b6*b6 + K::from_int(9)*b2*b4*b6
        }

        pub fn is_on_curve(&self, point : &ProjKPoint<K>) -> bool{
            use ProjKPoint::*;
            // y^2 + a_1 xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6
            match *point{
                InfPoint => true,
                FinPoint(x, y) =>{
                    y*y + self.a_1*x*y + self.a_3 * y == x*x*x + self.a_2*x*x + self.a_4*x + self.a_6
                }
            }
        }

        pub fn neg_point(&self, point : ProjKPoint<K>) -> ProjKPoint<K>{
            assert!(self.is_on_curve(&point));
            use ProjKPoint::*;
            match point {
                InfPoint => InfPoint,
                FinPoint(x, y) => FinPoint(x, -y - self.a_1*x - self.a_3),
            }
        }

        pub fn add_points(&self, point1 : ProjKPoint<K>, point2 : ProjKPoint<K>) -> ProjKPoint<K>{
            assert!(self.is_on_curve(&point1));
            assert!(self.is_on_curve(&point2));

            use ProjKPoint::*;

            let a_1 = self.a_1;
            let a_2 = self.a_2;
            let a_3 = self.a_3;
            let a_4 = self.a_4;
            let a_6 = self.a_6;

            if point1 == self.neg_point(point2){
                return InfPoint;
            }
            match (point1, point2){
                (InfPoint, _) => point2,
                (_, InfPoint) => point1,
                (FinPoint(x1, y1), FinPoint(x2, y2)) => {
                    let (lambda, nu) = 
                        if x1 != x2 {
                            ((y2-y1)/(x2-x1), (y1*x2 - y2*x1)/(x2-x1))
                        }else{
                            ((K::from_int(3)*x1*x1 + K::from_int(2)*a_2*x1 + a_4 - a_1*y1)/(K::from_int(2)*y1 + a_1*x1 + a_3), 
                             (-x1*x1*x1 + a_4*x1 + K::from_int(2)*a_6 - a_3*y1)/(K::from_int(2)*y1 + a_1*x1 + a_3))
                        };
                    let x3 = lambda*lambda +a_1*lambda - a_2 - x1 - x2;
                    let y3 = -(lambda + a_1)*x3 - nu - a_3;
                    FinPoint(x3, y3)
                }
            }
        }

        pub fn scalar_mult(&self, n : i32, point : ProjKPoint<K>) -> ProjKPoint<K>{
            assert!(self.is_on_curve(&point));
            if n == 0{
                return InfPoint;
            }

            if n < 0{
                return self.scalar_mult(-n, self.neg_point(point));
            }

            use ProjKPoint::*;
            let mut p1 = point;
            let mut p2 = self.add_points(p1, p1);

            let mut logm = -1; // -1 is here in order to ignore the first bit (included in p2 already)
            let mut m = n;
            while m != 0{
                m >>= 1;
                logm += 1;
            }

            while logm >= 1{
                let bit = (n&(1<<(logm-1)))>>(logm-1); // the current bit
                logm -= 1;
                if bit == 0{
                    p2 = self.add_points(p1, p2);
                    p1 = self.add_points(p1, p1);
                }else{
                    p1 = self.add_points(p1, p2);
                    p2 = self.add_points(p2, p2);
                }
            }
            p1
        }
}

#[cfg(test)]
mod test;