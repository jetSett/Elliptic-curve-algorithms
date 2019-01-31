use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign};
use std::marker::Sized;

pub type Integer = i128;

pub trait Field : Sized + 
                  Add<Output=Self> +
                  Sub<Output=Self> +
                  Mul<Output=Self> +
                  Div<Output=Self> +
                  Neg<Output=Self> +
                  AddAssign +
                  PartialEq +
                  Copy +
                  {
                    fn from_int(n : Integer) -> Self;
                  }
