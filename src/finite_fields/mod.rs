// p = 29

extern crate gmp;

use std::ops::{Add, Sub, Mul, Div};

use std::marker::PhantomData;

pub type Integer = i32;



pub trait IntegerAsType{
    fn value() -> Integer;
}

pub struct Fp<N : IntegerAsType>{
    repr : Integer,
    _phantom : PhantomData<N>,
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

        pub fn print(&self){
            println!("{}", self.repr);
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


impl<N> Add for &Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn add(self, other: &Fp<N>) -> Fp<N>{
            Fp::<N>::new((self.repr + other.repr) % N::value())
        }
}

impl<N> Add for Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn add(self, other: Fp<N>) -> Fp<N>{
            Fp::<N>::new((self.repr + other.repr) % N::value())
        }
}


impl<N> Sub for &Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn sub(self, other: &Fp<N>) -> Fp<N>{
            Fp::<N>::new((self.repr - other.repr + N::value()) % N::value())
        }
}

impl<N> Sub for Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn sub(self, other: Fp<N>) -> Fp<N>{
            Fp::<N>::new((self.repr - other.repr + N::value()) % N::value())
        }
}


impl<N> Mul for &Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn mul(self, other: &Fp<N>) -> Fp<N>{
            Fp::<N>::new((self.repr*other.repr) % N::value())
        }
}

impl<N> Mul for Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn mul(self, other: Fp<N>) -> Fp<N>{
            Fp::<N>::new((self.repr*other.repr) % N::value())
        }
}


impl<N> Div for &Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn div(self, other: &Fp<N>) -> Fp<N>{
            self*&(other.inv())
        }
}

impl<N> Div for Fp<N>
    where N : IntegerAsType{
        type Output = Fp<N>;

        fn div(self, other: Fp<N>) -> Fp<N>{
            self*(other.inv())
        }
}

impl<N> PartialEq for Fp<N>
    where N : IntegerAsType{
        fn eq(&self, other: &Fp<N>) -> bool{
            self.repr == other.repr
        }
}