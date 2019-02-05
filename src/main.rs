#![feature(type_alias_enum_variants)]

#[macro_use] pub mod finite_fields;
pub mod elliptic_curves;
pub mod field;

mod csidh;

use csidh::PublicKey;
use elliptic_curves::EllipticCurve;

fn main() {
    let mut pk = PublicKey::new(0);

    println!("Before class group action: {} - valid={}", pk, csidh::verify_public_key(pk));

    //println!("{}", EllipticCurve::new_montgomery(pk).j_invariant());

    let e = vec![1, 1, 1, 0, 0, 0];
    for _ in 1..10{
        let a = csidh::class_group_action(pk, e.clone());
        println!("{} - {}", a, EllipticCurve::new_montgomery(a).j_invariant());
        print!("\n");
    }

    // println!("{}", EllipticCurve::new_montgomery(pk).j_invariant());

    // println!("After class group action: {} - valid={}", pk, csidh::verify_public_key(pk));

    // let e1 = vec![1, -10, 0, -3, 1, -1];

    // pk = csidh::class_group_action(pk, e1);

    // println!("After class group action: {} - valid={}", pk, csidh::verify_public_key(pk));
}
