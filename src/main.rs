#![feature(type_alias_enum_variants)]

#[macro_use] pub mod finite_fields;
pub mod elliptic_curves;
pub mod field;

mod csidh;

use csidh::PublicKey;

fn main() {
    let mut pk = PublicKey::new(0);

    println!("Before class group action: valid={}", csidh::verify_public_key(pk));

    let e = vec![3, -5, 1, 2, -8, 3];

    pk = csidh::class_group_action(pk, e);

    println!("After class group action: valid={}", csidh::verify_public_key(pk));
}
