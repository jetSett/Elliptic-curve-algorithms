use crate::finite_fields::*;

use super::*;

use crate::field::IntegerTrait;

#[derive(Clone, Copy, Debug)]
pub struct UnsignedProjPoint<K : FiniteField>{
    pub x: K,
    pub z: K,
}

impl<K : FiniteField> UnsignedProjPoint<K> {
    pub fn infinite_point() -> UnsignedProjPoint<K>{
        UnsignedProjPoint{
            x: K::from_int(1), 
            z: K::from_int(0),
        }
    }

    pub fn finite_point(x : K) -> UnsignedProjPoint<K>{
        UnsignedProjPoint{
            x,
            z: K::from_int(1),
        }
    }

    pub fn order_two() -> UnsignedProjPoint<K>{
        UnsignedProjPoint{
            x: K::from_int(0),
            z: K::from_int(1),
        }
    }

    pub fn normalize(self) -> UnsignedProjPoint<K>{
        if self.z == K::from_int(0){
            Self::infinite_point()
        }else{
            Self::finite_point(self.x/self.z)
        }
    }
}

impl<K> PartialEq for UnsignedProjPoint<K> where K : FiniteField {
    fn eq(&self, other: &UnsignedProjPoint<K>) -> bool{
        match self.z == K::from_int(0){
            true => &other.z == &K::from_int(0),
            false => {
                if &other.z==&K::from_int(0){
                    false
                }else{
                    other.x.clone()*self.z.clone() == self.x.clone()*other.z.clone()
                }
            }
        }
    }
}

impl<K> fmt::Display for UnsignedProjPoint<K>
    where K : fmt::Display + FiniteField{
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "[{} : {}]", self.x, self.z)
        }
}



impl<K> EllipticCurve<K>
    where K : FiniteField + fmt::Display{
    
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
            let sample_element = | |{ K::new(K::Integer::sample_uniform(&K::Integer::from(0), &(K::cardinal()-K::Integer::from(1)))) };
            let mut x = sample_element();

            let f = |x : K| {
                x.clone()*x.clone()*x.clone() + self.a_2.clone()*x.clone()*x.clone() + self.a_4.clone()*x + self.a_6.clone()
            };
            while f(x.clone()).legendre_symbol() != 1{
                x = sample_element();
            }
            ProjKPoint::FinPoint(x.clone(), f(x).square_root())
        }

        pub fn sample_unsigned(&self) -> UnsignedProjPoint<K>{
            Self::unsigne_point(self.sample_point())
        }

        pub fn unsigne_point(p : ProjKPoint<K>) -> UnsignedProjPoint<K>{
            match p{
                ProjKPoint::InfPoint => UnsignedProjPoint::infinite_point(),
                ProjKPoint::FinPoint(x, _) => UnsignedProjPoint::finite_point(x)
            }
        }

        pub fn x_add(&self, p : UnsignedProjPoint<K>, q : UnsignedProjPoint<K>, p_minus_q : UnsignedProjPoint<K>) -> UnsignedProjPoint<K>{
            let u = (p.x.clone()-p.z.clone())*(q.x.clone()+q.z.clone());
            let v = (p.x+p.z)*(q.x-q.z);

            let x = u.clone()+v.clone();
            let z = u-v;

            UnsignedProjPoint{
                x: p_minus_q.z*x.clone()*x,
                z: p_minus_q.x*z.clone()*z,
            }
        }

        pub fn x_dbl(&self, p : UnsignedProjPoint<K>) -> UnsignedProjPoint<K>{
            let a = self.a_2.clone();
            let mut q = p.x.clone() + p.z.clone();
            q = q.clone()*q;

            let mut r = p.x - p.z;
            r = r.clone()*r;

            let s = q.clone() - r.clone();

            UnsignedProjPoint{
                x: q*r.clone(),
                z: s.clone()*(r + s*(a + K::from_int(2))/K::from_int(4) )
            }
        }

        pub fn scalar_mult_unsigned(&self, n : K::Integer, point : UnsignedProjPoint<K>) -> UnsignedProjPoint<K>{

            if n == K::Integer::from(0){
                return UnsignedProjPoint::infinite_point();
            }

            if n < K::Integer::from(0){
                return self.scalar_mult_unsigned(-n, point);
            }
            
            let mut logm = 0;
            let mut m = n.clone();
            while m != K::Integer::from(0){
                m >>= 1;
                logm += 1;
            }

            let mut x0 = UnsignedProjPoint::infinite_point();
            let mut x1 = point.clone();

            let cond_swap = |b: K, x0: UnsignedProjPoint<K>, x1: UnsignedProjPoint<K>| {
                (UnsignedProjPoint{
                    x: (K::from_int(1)-b.clone())*x0.x.clone() + b.clone()*x1.x.clone(), 
                    z: (K::from_int(1)-b.clone())*x0.z.clone() + b.clone()*x1.z.clone(), 
                }, 
                UnsignedProjPoint{
                    x: (K::from_int(1)-b.clone())*x1.x.clone() + b.clone()*x0.x.clone(),
                    z: (K::from_int(1)-b.clone())*x1.z.clone() + b.clone()*x0.z.clone(),
                })
            };

            while logm >= 1{
                let bit = K::new((n.clone()&(K::Integer::from(1)<<(logm-1)))>>(logm-1)); // the current bit
                logm -= 1;

                let  (mut a, mut b) = cond_swap(bit.clone(), x0, x1);

                b = self.x_add(a.clone(), b, point.clone());
                a = self.x_dbl(a);

                let (_a, _b) = cond_swap(bit, a, b);
                x0 = _a;
                x1 = _b;
            }
            x0            
        }


}

#[cfg(test)]
mod test;