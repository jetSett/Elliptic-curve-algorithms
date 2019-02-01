use rand::Rng;

use crate::finite_fields::*;

use super::*;

pub enum ProjFpPoint<K : FiniteField>{
    InfPoint,
    FinPoint(K, bool)
}
#[derive(Clone, Copy)]
struct InternalProjPoint<K : FiniteField>{
    x: K,
    z: K,
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

        pub fn decompress_point(&self, p : ProjFpPoint<K>) -> ProjKPoint<K>{
            if !self.left_side_empty(){
                panic!("Point compression must be used with only y^2");
            }
            match p{
                ProjFpPoint::InfPoint => ProjKPoint::InfPoint,
                ProjFpPoint::FinPoint(x, y_sign) => {
                    let y = (x*x*x + self.a_2*x*x + self.a_4*x + self.a_6).square_root()* 
                            if y_sign {K::from_int(1)} else {K::from_int(-1)};
                    ProjKPoint::FinPoint(x, y)
                }
            }        
        }

        pub fn compress_point(&self, p : ProjKPoint<K>) -> ProjFpPoint<K>{
            match p{
                ProjKPoint::InfPoint => ProjFpPoint::InfPoint,
                ProjKPoint::FinPoint(x, y) => ProjFpPoint::FinPoint(x, y.sign())
            }
        }

        pub fn neg_compressed_point(&self, p : ProjFpPoint<K>) -> ProjFpPoint<K>{
            match p{
                ProjFpPoint::InfPoint => ProjFpPoint::InfPoint,
                ProjFpPoint::FinPoint(x, y) => ProjFpPoint::FinPoint(x, !y)
            }
        }

        pub fn add_compressed_point(&self, p1 : ProjFpPoint<K>, p2 : ProjFpPoint<K>) -> ProjFpPoint<K>{
            assert!(self.is_montgomery());
            let p1d = self.decompress_point(p1);
            let p2d = self.decompress_point(p2);

            let p3d = self.add_points(p1d, p2d);

            self.compress_point(p3d)
        }

        fn compressed_to_internal(p : ProjFpPoint<K>) -> InternalProjPoint<K>{
            match p{
                ProjFpPoint::InfPoint => InternalProjPoint{
                    x: K::from_int(1), z: K::from_int(0)
                },
                ProjFpPoint::FinPoint(x, _) => InternalProjPoint{
                    x, z: K::from_int(1)
                },
            }
        }

        fn internal_to_compressed(p : InternalProjPoint<K>, sign: bool) -> ProjFpPoint<K>{
            if p.z == K::from_int(0){
                ProjFpPoint::InfPoint
            }else{
                ProjFpPoint::FinPoint(p.x/p.z, sign)
            }
        }


        fn x_add(&self, p : InternalProjPoint<K>, q : InternalProjPoint<K>, p_minus_q : InternalProjPoint<K>) -> InternalProjPoint<K>{
            let u = (p.x-p.z)*(q.x+q.z);
            let v = (p.x+p.z)*(q.x-q.z);

            InternalProjPoint{
                x: p_minus_q.z*(u+v)*(u+v),
                z: p_minus_q.x*(u-v)*(u-v),
            }
        }

        fn x_dbl(&self, p : InternalProjPoint<K>) -> InternalProjPoint<K>{
            let a = self.a_4;
            let q = (p.x + p.z)*(p.x + p.z);
            let r = (p.x - p.z)*(p.x - p.z);
            let s = q-r;
            InternalProjPoint{
                x: q*r,
                z: s*(r + s*(a + K::from_int(2))/K::from_int(4))
            }
        }

        pub fn scalar_mult_compressed(&self, n : Integer, point : ProjFpPoint<K>) -> ProjFpPoint<K>{
            use ProjFpPoint::*;

            if n == 0{
                return InfPoint;
            }

            if n < 0{
                return self.scalar_mult_compressed(-n, self.neg_compressed_point(point));
            }
            

            if let FinPoint(_, sign_y) = point{
                let mut logm = 0;
                let mut m = n;
                while m != 0{
                    m >>= 1;
                    logm += 1;
                }

                let internal_p = Self::compressed_to_internal(point);

                let (mut x0, mut x1) = (Self::compressed_to_internal(ProjFpPoint::InfPoint), internal_p);

                let cond_swap = |b: K, x0: InternalProjPoint<K>, x1: InternalProjPoint<K>| {
                    (InternalProjPoint{
                        x: (K::from_int(1)-b)*x0.x + b*x1.x, 
                        z: (K::from_int(1)-b)*x0.z + b*x1.z, 
                    }, 
                    InternalProjPoint{
                        x: (K::from_int(1)-b)*x1.x + b*x0.x,
                        z: (K::from_int(1)-b)*x1.z + b*x0.z,
                    })
                };

                while logm >= 1{
                    let bit = K::from_int((n&(1<<(logm-1)))>>(logm-1)); // the current bit
                    logm -= 1;

                    let  (mut a, mut b) = cond_swap(bit, x0, x1);
                    a = self.x_dbl(a);
                    b = self.x_add(a, b, internal_p);
                    let (_a, _b) = cond_swap(bit, a, b);
                    x0 = a;
                    x1 = b;
                }
                
                Self::internal_to_compressed(x0, sign_y)
            }else{
                InfPoint
            }

        }


}

#[cfg(test)]
mod test;