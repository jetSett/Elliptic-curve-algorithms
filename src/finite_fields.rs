use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, MulAssign, DivAssign, SubAssign};
use std::clone::Clone;
use std::fmt;

use std::marker::PhantomData;

use rand::Rng;

use crate::field::{Field, Integer};

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
            use super::*;

            #[derive(Debug)]
            pub struct TypeInt{}

            impl IntegerAsType for TypeInt{
                fn value() -> Integer{
                    $p
                }
            }
        }


        pub type $name = Fp<$m::TypeInt>;
    }
}

pub trait FiniteField : Field{
    fn legendre_symbol(&self) -> i8;

    fn exp(a : Self, n : Integer) -> Self;

    fn square_root(self) -> Self;

    fn cardinal() -> Integer;

    fn sign(&self) -> bool; // true = +, false = -, + = closest to 0
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

impl<N> AddAssign for Fp<N>
    where N : IntegerAsType{
        fn add_assign(&mut self, other: Fp<N>){
            *self = *self + other;
        }
}

impl<N> MulAssign for Fp<N>
    where N : IntegerAsType{
        fn mul_assign(&mut self, other: Fp<N>){
            *self = *self*other;
        }
}

impl<N> DivAssign for Fp<N>
    where N : IntegerAsType{
        fn div_assign(&mut self, other: Fp<N>){
            *self = *self/other;
        }
}

impl<N> SubAssign for Fp<N>
    where N : IntegerAsType{
        fn sub_assign(&mut self, other: Fp<N>){
            *self = *self-other;
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

impl<N> Field for Fp<N>
    where N : IntegerAsType {
        fn from_int(n : Integer) -> Fp<N> {
            Fp::<N>::new(n)
        }
    }

impl<N> FiniteField for Fp<N>
    where N : IntegerAsType{
        fn legendre_symbol(&self) -> i8{
            let mut a = self.repr;
            let mut m = N::value();
            let mut t : i8= 1;
            while a != 0{
                while a%2 == 0{
                    a = a/2;
                    if m%8 == 3 || m%8 == 5{
                        t = -t;
                    }
                }
                let _c = a;
                a = m;
                m = _c;

                if a%4 == 3 && m%4 == 3{
                    t = -t;
                }
                a = a%m;
            }
            if m==1{
                t
            }else{
                0
            }
        }

        fn cardinal() -> Integer{
            N::value()
        }

        fn exp(a : Fp<N>, n : Integer) -> Fp<N>{
            if n == 0{
                return Fp::<N>::new(1);
            }
            if n<0{
                return Fp::<N>::exp(a.inv(), -n);
            }
            if n%2 == 0{
                Fp::<N>::exp(a*a, n/2)
            }else{
                a*Fp::<N>::exp(a*a, (n-1)/2)
            }
        }

        fn square_root(self) -> Fp<N>{

            let exp = Fp::exp;

            let p = N::value();
            let a = self;

            if p%8 == 3 || p%8 == 5 || p%8 == 7{
                return exp(a, (p+1)/4);
            }

            if p%8 == 5{
                let mut x = exp(a, (p+3)/8);
                let c = x*x;
                if c != self{
                    x = x*exp(Fp::<N>::new(2), (p-1)/4);
                }
                return x;
            }

            // Here p = 1 mod 8
            let mut rng = rand::thread_rng();
            let mut d = Fp::<N>::new(rng.gen_range(2, p-1));

            while d.legendre_symbol() != -1{
                d = Fp::<N>::new(rng.gen_range(2, p-1));
            }

            let mut t = p-1;
            let mut s = 0;
            while t%2 == 0{ // represent p-1 = t*2^s
                s += 1;
                t >>= 1;
            }

            let big_a = exp(a, t);
            let big_d = exp(d, t);

            let mut m = 0;

            let mut exponent = 1 << s-1;

            for i in 0..s{
                if exp(big_a*exp(big_d, m), exponent) == Fp::new(-1){
                    m += 1<<i;
                }
                exponent >>= 1
            }
            
            let sqr = exp(a, (t+1)/2)*exp(big_d, m/2);
            if sqr.sign(){
                sqr
            }else{
                -sqr
            }
        }

        fn sign(&self) -> bool {
            self.repr <= (N::value()-1)/2 // true if self is closer to 0 (0 is positive)
        }
}



#[cfg(test)]
mod test;