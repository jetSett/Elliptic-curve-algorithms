#![feature(type_alias_enum_variants)]

use time;

#[macro_use] pub mod finite_fields;
pub mod elliptic_curves;
pub mod field;

mod csidh;

use csidh::*;
use elliptic_curves::EllipticCurve;

fn main() {
    let gen_param = 10;

    let mut current : u64;

    current = time::precise_time_ns();
    let (pk_a, sk_a) = sample_keys(gen_param);
    let time_gen_a = (time::precise_time_ns()-current)/1000000;
    println!("Alice's public key: {} ({} ms)", pk_a, time_gen_a);

    current = time::precise_time_ns();
    let (pk_b, sk_b) = sample_keys(gen_param);
    let time_gen_b = (time::precise_time_ns()-current)/1000000;
    println!("Bob's public key: {} ({} ms)", pk_b, time_gen_b);

    current = time::precise_time_ns();
    let shared_a = class_group_action(pk_b, sk_a);
    let time_shared_a = (time::precise_time_ns()-current)/1000000;
    println!("Alice's shared secret: {} ({} ms)", shared_a, time_shared_a);


    current = time::precise_time_ns();
    let shared_b = class_group_action(pk_a, sk_b);
    let time_shared_b = (time::precise_time_ns()-current)/1000000;
    println!("Bob's shared secret: {} ({} ms)", shared_b, time_shared_b);

}
