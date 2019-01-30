use std::ops::{Add, Sub, Mul, Div, Neg};
use std::marker::Sized;

pub type Integer = i64;

pub trait FieldValues<K>{
    fn from_int(n : Integer) -> K;
}

pub trait Field : Sized + 
                  Add<Output=Self> +
                  Sub<Output=Self> +
                  Mul<Output=Self> +
                  Div<Output=Self> +
                  Neg<Output=Self> +
                  PartialEq +
                  Copy +
                  {}
