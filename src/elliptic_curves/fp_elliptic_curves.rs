use rand::Rng;

use crate::finite_fields::*;

use super::*;

impl<N> EllipticCurve<Fp<N>> where N : IntegerAsType {

    pub fn efficient_sample_point(&self) -> ProjKPoint<Fp<N>>{
        if self.a_1 != Fp::new(0) || self.a_3 != Fp::new(0){
            panic!("efficient_sample_point point must be used with curves with only y^2");
        }

        let mut rng = rand::thread_rng();
        let mut x = Fp::new(rng.gen_range(0, N::value()-1));

        let f = |x| {
            x*x*x + self.a_2*x*x + self.a_4*x + self.a_6
        };
        while f(x).legendre_symbol() != 1{
            x = Fp::new(rng.gen_range(0, N::value()-1));
        }
        ProjKPoint::FinPoint(x, f(x).square_root())
    }


}
