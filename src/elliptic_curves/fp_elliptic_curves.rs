use rand::Rng;

use crate::finite_fields::*;

use super::*;

pub enum ProjFpPoint<K : FiniteField>{
    InfPoint,
    FinPoint(K, bool)
}

impl<K> EllipticCurve<K>
    where K : FiniteField{
    
        fn left_side_empty(&self) -> bool{
            self.a_1 == K::from_int(0) && self.a_3 == K::from_int(0)
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
            let p1d = self.decompress_point(p1);
            let p2d = self.decompress_point(p2);

            let p3d = self.add_points(p1d, p2d);

            self.compress_point(p3d)
        }

}
