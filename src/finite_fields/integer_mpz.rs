use std::ops::{Add, Sub, Mul};

use gmp::mpz::Mpz;

use super::*;

use crate::field::{Field};

impl<N : IntegerAsType<Mpz>> FiniteField for Fp<N, Mpz>{
    fn exp(a : Fp<N, Mpz>, n : Mpz) -> Fp<N, Mpz>{
        Fp::new(a.repr.powm(&n, &N::value()))
    }
}

impl<N : IntegerAsType<Mpz>> num_traits::ops::inv::Inv for Fp<N, Mpz>{
    fn inv(self) -> Fp<N, Mpz>{
        assert!(self.repr != Mpz::from(0));
        Fp::new(self.repr.invert(&N::value()).unwrap())
    }
}

impl<N : IntegerAsType<Mpz>> Add for &Fp<N, Mpz>{
    fn add(self, other: &Fp<N, Mpz>) -> Fp<N, Mpz>{
        Fp::new((&self.repr+&other.repr)%N::value())
    }
}

impl<N : IntegerAsType<Mpz>> Sub for &Fp<N, Mpz>{
    fn sub(self, other: &Fp<N, Mpz>) -> Fp<N, Mpz>{
        if self.repr > other.repr{
            Fp::new(&self.repr - &other.repr)
        }else{
            Fp::new(&self.repr - &other.repr + N::value())
        }
    }
}

impl<N : IntegerAsType<Mpz>> Mul for &Fp<N, Mpz>{
    fn mul(self, other: &Fp<N, Mpz>) -> Fp<N, Mpz>{
        Fp::new((&self.repr*&other.repr)%N::value())
    }
}


