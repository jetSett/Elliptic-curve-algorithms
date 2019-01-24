use std::ops::{Add, Sub, Mul, Div, Neg};
use std::clone::Clone;
use std::fmt;

use std::marker::PhantomData;

pub type Integer = i64;

pub trait IntegerAsType{
    fn value() -> Integer;
}

#[derive(Debug)]
pub struct Fp<N : IntegerAsType>{
    repr : Integer,
    _phantom : PhantomData<N>,
}

#[macro_export]
macro_rules! declare_finite_field {
    ($name: ident, $p: expr, $m:ident) => {
        mod $m{
            #[derive(Debug)]
            pub struct TypeInt{}
        }

        impl IntegerAsType for $m::TypeInt{
            fn value() -> Integer{
                $p
            }
        }

        pub type $name = Fp<$m::TypeInt>;
    }
}

impl<N> fmt::Display for Fp<N>
    where N : IntegerAsType{
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{}", self.repr)
        }
}

impl<N> Fp<N>
    where N : IntegerAsType{
        pub fn new(x : Integer) -> Fp<N>{
            Fp::<N>{
                repr: (x%N::value() + N::value())%N::value(),
                _phantom: PhantomData,
            }
        }

        pub fn zero() -> Fp<N>{
            Fp::new(0)
        }

        pub fn inv(&self) -> Fp<N>{
            assert!(self.repr != 0);
            let (mut t, mut new_t) = (0, 1);
            let (mut r, mut new_r) = (N::value(), self.repr);
            
            while new_r != 0{
                let quotient = r / new_r;
                
                let old_t = t;
                t = new_t;
                new_t = old_t - quotient*new_t;

                let old_r = r;
                r = new_r;
                new_r = old_r - quotient*new_r;
            }
            if t < 0{
                t += N::value();
            }
            Fp::<N>::new(t)
        }
}

impl<N> Clone for Fp<N>
    where N : IntegerAsType{
        fn clone(&self) -> Fp<N>{
            Fp::<N>::new(self.repr)
        }
}

impl<N> Copy for Fp<N> where N : IntegerAsType{}

impl<N> Add for Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn add(self, other: Fp<N>) -> Fp<N>{
            Fp::<N>::new((self.repr + other.repr) % N::value())
        }
}

impl<N> Sub for Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn sub(self, other: Fp<N>) -> Fp<N>{
            Fp::<N>::new((self.repr - other.repr + N::value()) % N::value())
        }
}

impl<N> Mul for Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn mul(self, other: Fp<N>) -> Fp<N>{
            Fp::<N>::new((self.repr*other.repr) % N::value())
        }
}

impl<N> Div for Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn div(self, other: Fp<N>) -> Fp<N>{
            self*(other.inv())
        }
}

impl<N> Neg for Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn neg(self) -> Fp<N>{
            Fp::<N>::new(-self.repr)
        }
}

impl<N> PartialEq for Fp<N>
    where N : IntegerAsType{
        fn eq(&self, other: &Fp<N>) -> bool{
            self.repr == other.repr
        }
}

mod test;
