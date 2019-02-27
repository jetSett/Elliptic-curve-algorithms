use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, MulAssign, DivAssign, SubAssign};
use std::clone::Clone;
use std::fmt;

use std::marker::PhantomData;

use crate::field::{Field, IntegerTrait};

pub trait IntegerAsType<Integer : IntegerTrait>{
    fn value() -> Integer;
}

#[derive(Debug, Hash)]
pub struct Fp<N : IntegerAsType<Integer>, Integer : IntegerTrait>{
    repr : Integer,
    _phantom : PhantomData<N>,
}

#[macro_export]
macro_rules! declare_finite_field {
    ($name: ident, $integer: ident, $p: expr, $m:ident) => {
        mod $m{
            use super::*;

            #[derive(Debug, Hash)]
            pub struct TypeInt{}

            impl IntegerAsType<$integer> for TypeInt{
                fn value() -> $integer{
                    $p
                }
            }
        }


        pub type $name = Fp<$m::TypeInt, $integer>;
    }
}

pub trait FiniteField : Field{

    fn legendre_symbol(&self) -> i8;

    fn exp(a : Self, n : Self::Integer) -> Self;

    fn square_root(self) -> Self;

    fn cardinal() -> Self::Integer;

    fn sign(&self) -> bool; // true = +, false = -, + = closest to 0
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> fmt::Display for Fp<N, Integer>{
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{}", self.repr)
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Fp<N, Integer>{

        pub fn inv(&self) -> Fp<N, Integer>{
            assert!(self.repr != Integer::from(0));
            let (mut t, mut new_t) = (Integer::from(0), Integer::from(1));
            let (mut r, mut new_r) = (N::value(), self.repr.clone());
            
            while new_r != Integer::from(0){
                let quotient = r.clone() / new_r.clone();
                
                let old_t = t;
                t = new_t.clone();
                new_t = old_t - quotient.clone()*new_t;

                let old_r = r;
                r = new_r.clone();
                new_r = old_r - quotient*new_r;
            }
            if t < Integer::from(0){
                t += N::value();
            }
            Fp::new(t)
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Clone for Fp<N, Integer>{
        fn clone(&self) -> Fp<N, Integer>{
            Fp::new(self.repr.clone())
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Add for Fp<N, Integer>{
        type Output = Fp<N, Integer>;

        fn add(self, other: Fp<N, Integer>) -> Fp<N, Integer>{
            Fp::new((self.repr + other.repr) % N::value())
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Add for &Fp<N, Integer>{
        type Output = Fp<N, Integer>;

        fn add(self, other: &Fp<N, Integer>) -> Fp<N, Integer>{
            Fp::new((self.repr.clone() + other.repr.clone()) % N::value())
        }
}


impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> AddAssign for Fp<N, Integer>{
        fn add_assign(&mut self, other: Fp<N, Integer>){
            *self = self.clone() + other;
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> MulAssign for Fp<N, Integer>{
        fn mul_assign(&mut self, other: Fp<N, Integer>){
            *self = self.clone()*other;
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> DivAssign for Fp<N, Integer>{
        fn div_assign(&mut self, other: Fp<N, Integer>){
            *self = self.clone()/other;
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> SubAssign for Fp<N, Integer>{
        fn sub_assign(&mut self, other: Fp<N, Integer>){
            *self = self.clone()-other;
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Sub for Fp<N, Integer>{
        type Output = Fp<N, Integer>;

        fn sub(self, other: Fp<N, Integer>) -> Fp<N, Integer>{
            Fp::new((self.repr - other.repr + N::value()) % N::value())
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Sub for &Fp<N, Integer>{
        type Output = Fp<N, Integer>;

        fn sub(self, other: &Fp<N, Integer>) -> Fp<N, Integer>{
            Fp::new((self.repr.clone() - other.repr.clone() + N::value()) % N::value())
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Mul for Fp<N, Integer>{
        type Output = Fp<N, Integer>;

        fn mul(self, other: Fp<N, Integer>) -> Fp<N, Integer>{
            Fp::new((self.repr*other.repr) % N::value())
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Mul for &Fp<N, Integer>{
        type Output = Fp<N, Integer>;

        fn mul(self, other: &Fp<N, Integer>) -> Fp<N, Integer>{
            Fp::new((self.repr.clone()*other.repr.clone()) % N::value())
        }
}


impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Div for Fp<N, Integer>{
        type Output = Fp<N, Integer>;

        fn div(self, other: Fp<N, Integer>) -> Fp<N, Integer>{
            self*(other.inv())
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Neg for Fp<N, Integer>{
        type Output = Fp<N, Integer>;

        fn neg(self) -> Fp<N, Integer>{
            Fp::new(-self.repr)
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> PartialEq for Fp<N, Integer>{
        fn eq(&self, other: &Fp<N, Integer>) -> bool{
            self.repr == other.repr
        }
}

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> Field for Fp<N, Integer> {
        
        type Integer=Integer;

        fn new(x : Integer) -> Fp<N, Integer>{
            Fp{
                repr: (x%N::value() + N::value())%N::value(),
                _phantom: PhantomData,
            }
        }

        fn from_int(n : i32) -> Fp<N, Integer> {
            Fp::new(Integer::from(n))
        }
    }

impl<N : IntegerAsType<Integer>, Integer : IntegerTrait> FiniteField for Fp<N, Integer>{

    fn legendre_symbol(&self) -> i8{
        let mut a = self.repr.clone();
        let mut m = N::value();
        let mut t : i8 = 1;
        while a != Integer::from(0){
            while a.clone()%2 == Integer::from(0){
                a = a/Integer::from(2);
                if m.clone()%8 == Integer::from(3) || m.clone()%8 == Integer::from(5){
                    t = -t;
                }
            }
            let _c = a;
            a = m;
            m = _c;

            if a.clone()%4 == Integer::from(3) && m.clone()%4 == Integer::from(3){
                t = -t;
            }
            a = a%m.clone();
        }
        if m==Integer::from(1){
            t
        }else{
            0
        }
    }

    fn cardinal() -> Integer{
        N::value()
    }

    fn exp(a : Fp<N, Integer>, n : Integer) -> Fp<N, Integer>{
        if n == Integer::from(0){
            return Fp::from_int(1);
        }
        if n<Integer::from(0){
            return Fp::exp(a.inv(), -n);
        }
        if n.clone()%2 == Integer::from(0){
            Fp::exp(&a*&a, n/Integer::from(2))
        }else{
            &a*&Fp::exp(&a*&a, (n-Integer::from(1))/Integer::from(2))
        }
    }

    fn square_root(self) -> Fp<N, Integer>{

        let exp = Fp::exp;

        let p = N::value();
        let a = self.clone();

        if p.clone()%8 == Integer::from(3) || p.clone()%8 == Integer::from(7){
            return exp(a, (p.clone()+Integer::from(1))/Integer::from(4));
        }

        if p.clone()%8 == Integer::from(5){
            let mut x = exp(a, (p.clone()+Integer::from(3))/Integer::from(8));
            let c = &x*&x;
            if c != self{
                x = &x*&exp(Fp::from_int(2), (p.clone()-Integer::from(1))/Integer::from(4));
            }
            return x;
        }

        // Here p = 1 mod 8
        let mut d = Fp::new(Integer::sample_uniform(&Integer::from(2), &(p.clone()-Integer::from(1))));

        while d.legendre_symbol() != -1{
            d = Fp::new(Integer::sample_uniform(&Integer::from(2), &(p.clone()-Integer::from(1))));
        }

        let mut t = p-Integer::from(1);
        let mut s = 0;
        while t.clone()%2 == Integer::from(0){ // represent p-1 = t*2^s
            s += 1;
            t >>= 1;
        }

        let big_a = exp(a.clone(), t.clone());
        let big_d = exp(d, t.clone());

        let mut m = Integer::from(0);

        let mut exponent : Integer = Integer::from(1) << (s-1);

        for i in 0..s{
            if exp(big_a.clone()*exp(big_d.clone(), m.clone()), exponent.clone()) == Fp::from_int(-1){
                m += Integer::from(1)<<i;
            }
            exponent >>= 1
        }
        
        let sqr = exp(a, (t+Integer::from(1))/Integer::from(2))*exp(big_d, m/Integer::from(2));
        if sqr.sign(){
            sqr
        }else{
            -sqr
        }
    }

    fn sign(&self) -> bool {
        self.repr <= (N::value()-Integer::from(1))/Integer::from(2) // true if self is closer to 0 (0 is positive)
    }
}



#[cfg(test)]
mod test;