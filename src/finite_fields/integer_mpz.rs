use gmp::mpz::Mpz;

use crate::field::{IntegerTrait, Field};
use super::{Fp, FiniteField, IntegerAsType};

type Integer = Mpz;

impl<N : IntegerAsType<Mpz>> FiniteField for Fp<N, Mpz>{
    fn exp(a : Fp<N, Integer>, n : Integer) -> Fp<N, Integer>{
        Fp::new(a.repr.powm(&n, &N::value()))
    }
}

impl<N : IntegerAsType<Integer>> num_traits::ops::inv::Inv for Fp<N, Mpz>{
    fn inv(self) -> Fp<N, Integer>{
        assert!(self.repr != Integer::from(0));
        Fp::new(self.repr.invert(&N::value()).unwrap())
    }
}
