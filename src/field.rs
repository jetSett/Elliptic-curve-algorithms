use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, MulAssign, DivAssign, SubAssign};
use std::fmt::Display;
use std::marker::Sized;

pub type Integer = i128;

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
                  Copy +
                  Display
                  {
                    fn from_int(n : Integer) -> Self;
                  }
