#![feature(type_alias_enum_variants)]

#[macro_use] pub mod finite_fields;
pub mod elliptic_curves;
pub mod field;

mod csidh;

use finite_fields::*;
use elliptic_curves::*;

use field::*;

const P : Integer = 78539;

use csidh::PublicKey;

fn main() {
}
