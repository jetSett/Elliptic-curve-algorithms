use std::ops::{Add, Sub, Mul, 
              Div, Neg, AddAssign, 
              MulAssign, DivAssign, 
              SubAssign, Rem,
              Shr, ShrAssign,
              Shl, ShlAssign,
              BitAnd, BitOr};

use std::cmp::{PartialOrd};

use std::fmt::Display;
use std::marker::Sized;
use gmp::mpz::{Mpz};

pub trait Field : Sized + 
                  Add<Output=Self> +
                  Sub<Output=Self> +
                  Mul<Output=Self> +
                  Div<Output=Self> +
                  Neg<Output=Self> +
                  AddAssign +
                  MulAssign +
                  DivAssign +
                  SubAssign +
                  PartialEq +
                  Display +

                  Clone +

                  {
                    type Integer : IntegerTrait;
                    fn from_int(n : i32) -> Self;
                    fn new(n : Self::Integer) -> Self;
                  }

pub trait IntegerTrait : Sized + 
                  Add<Output=Self> +
                  Sub<Output=Self> +
                  Mul<Output=Self> +
                  Div<Output=Self> +
                  Neg<Output=Self> +
                  Rem<Self, Output=Self> +
                  Rem<u64, Output=Self> +
                  Shr<usize, Output=Self> +
                  Shl<usize, Output=Self> +
                  BitAnd<Output=Self> +
                  BitOr<Output=Self> +

                  AddAssign +
                  MulAssign +
                  DivAssign +
                  SubAssign +
                  ShrAssign<usize> +
                  ShlAssign<usize> +

                  PartialEq +
                  PartialOrd +

                  Clone +

                  Display + 
                  From<u32> + From<i32> + From<u32> + From<u64>{

                    fn sample_uniform(min : &Self, max : &Self) -> Self;

                  }


impl IntegerTrait for Mpz{
  fn sample_uniform(min : &Mpz, max : &Mpz) -> Mpz{
    let mut rng = gmp::rand::RandState::new();
    rng.urandom(&(max-min)) + min
  }
}
