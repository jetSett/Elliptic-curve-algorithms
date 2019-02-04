use rand::Rng;

use crate::finite_fields::*;

use super::*;

#[derive(Clone, Copy, Debug)]
pub struct UnsignedProjPoint<K : FiniteField>{
    x: K,
    z: K,
}

impl<K : FiniteField> UnsignedProjPoint<K> {
    fn infinite_point() -> UnsignedProjPoint<K>{
        UnsignedProjPoint{
            x: K::from_int(1), 
            z: K::from_int(0),
        }
    }

    fn finite_point(x : K) -> UnsignedProjPoint<K>{
        UnsignedProjPoint{
            x,
            z: K::from_int(1),
        }
    }
}

impl<K> PartialEq for UnsignedProjPoint<K> where K : FiniteField {
    fn eq(&self, other: &UnsignedProjPoint<K>) -> bool{
        match self.z == K::from_int(0){
            true => other.z == K::from_int(0),
            false => {
                if other.z==K::from_int(0){
                    false
                }else{
                    other.x/other.z == self.x/self.z
                }
            }
        }
    }
}

impl<K> EllipticCurve<K>
    where K : FiniteField{
    
        fn left_side_empty(&self) -> bool{
            self.a_1 == K::from_int(0) && self.a_3 == K::from_int(0)
        }

        pub fn is_montgomery(&self) -> bool{
            self.left_side_empty() && self.a_4 == K::from_int(1) && self.a_6 == K::from_int(0)
        }

        pub fn new_montgomery(a : K) -> EllipticCurve<K>{
            EllipticCurve::<K>{
                a_1: K::from_int(0),
                a_3: K::from_int(0),

                a_2: a,
                a_4: K::from_int(1),
                a_6: K::from_int(0),
            }
        }

        pub fn sample_point(&self) -> ProjKPoint<K>{
            if !self.left_side_empty(){
                panic!("sample_point point must be used with curves with only y^2");
            }

            let mut rng = rand::thread_rng();
            let mut x = K::from_int(rng.gen_range(0, K::cardinal()-1));

            let f = |x| {
                x*x*x + self.a_2*x*x + self.a_4*x + self.a_6
            };
            while f(x).legendre_symbol() != 1{
                x = K::from_int(rng.gen_range(0, K::cardinal()-1));
            }
            ProjKPoint::FinPoint(x, f(x).square_root())
        }

        pub fn unsigne_point(p : ProjKPoint<K>) -> UnsignedProjPoint<K>{
            match p{
                ProjKPoint::InfPoint => UnsignedProjPoint::infinite_point(),
                ProjKPoint::FinPoint(x, _) => UnsignedProjPoint::finite_point(x)
            }
        }

        fn x_add(&self, p : UnsignedProjPoint<K>, q : UnsignedProjPoint<K>, p_minus_q : UnsignedProjPoint<K>) -> UnsignedProjPoint<K>{
            let u = (p.x-p.z)*(q.x+q.z);
            let v = (p.x+p.z)*(q.x-q.z);

            UnsignedProjPoint{
                x: p_minus_q.z*(u+v)*(u+v),
                z: p_minus_q.x*(u-v)*(u-v),
            }
        }

        fn x_dbl(&self, p : UnsignedProjPoint<K>) -> UnsignedProjPoint<K>{
            let a = self.a_2;
            let q = (p.x + p.z)*(p.x + p.z);
            let r = (p.x - p.z)*(p.x - p.z);
            let s = q - r;
            UnsignedProjPoint{
                x: q*r,
                z: s*(r + s*(a + K::from_int(2))/K::from_int(4) )
            }
        }

        pub fn scalar_mult_unsigned(&self, n : Integer, point : UnsignedProjPoint<K>) -> UnsignedProjPoint<K>{

            if n == 0{
                return UnsignedProjPoint::infinite_point();
            }

            if n < 0{
                return self.scalar_mult_unsigned(-n, point);
            }
            
            let mut logm = 0;
            let mut m = n;
            // println!("{:b}", m);
            while m != 0{
                m >>= 1;
                logm += 1;
            }

            let mut x0 = UnsignedProjPoint::infinite_point();
            let mut x1 = point;

            let cond_swap = |b: K, x0: UnsignedProjPoint<K>, x1: UnsignedProjPoint<K>| {
                (UnsignedProjPoint{
                    x: (K::from_int(1)-b)*x0.x + b*x1.x, 
                    z: (K::from_int(1)-b)*x0.z + b*x1.z, 
                }, 
                UnsignedProjPoint{
                    x: (K::from_int(1)-b)*x1.x + b*x0.x,
                    z: (K::from_int(1)-b)*x1.z + b*x0.z,
                })
            };

            while logm >= 1{
                let bit = K::from_int((n&(1<<(logm-1)))>>(logm-1)); // the current bit
                // println!("{}", bit);
                logm -= 1;

                let  (mut a, mut b) = cond_swap(bit, x0, x1);
                a = self.x_dbl(a);
                b = self.x_add(a, b, point);
                let (_a, _b) = cond_swap(bit, a, b);
                x0 = a;
                x1 = b;
            }
            x0            
        }


}

#[cfg(test)]
mod test;