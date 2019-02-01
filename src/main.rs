#![feature(type_alias_enum_variants)]

#[macro_use] pub mod finite_fields;
pub mod elliptic_curves;
pub mod field;

mod csidh;

use finite_fields::*;
use elliptic_curves::*;

use field::*;

use csidh::PublicKey;

fn main() {
    let ell0 = EllipticCurve::new_montgomery(PublicKey::new(0));

    let pk = PublicKey::new(0);

    println!("{}", csidh::verify_public_key(pk));

    let e = vec![3, -5, 1, 2, -8, 3];

    csidh::class_group_action(pk, e);

}
