#![feature(type_alias_enum_variants)]

use time;

#[macro_use] pub mod finite_fields;
pub mod elliptic_curves;
pub mod field;

mod csidh;

use csidh::*;

fn main() {
    let gen_param = 10;

    let mut current : u64;

    let inst = CSIDHInstance{ 
        p: Integer::from(1021019),
        l: [Integer::from(3), Integer::from(5), Integer::from(7), Integer::from(11), Integer::from(13), Integer::from(17)]
    };

    current = time::precise_time_ns();
    let (pk_a, sk_a) = sample_keys(&inst, gen_param);
    let time_gen_a = (time::precise_time_ns()-current)/1000000;
    println!("Alice's public key: {} ({} ms)", pk_a, time_gen_a);

    current = time::precise_time_ns();
    let (pk_b, sk_b) = sample_keys(&inst, gen_param);
    let time_gen_b = (time::precise_time_ns()-current)/1000000;
    println!("Bob's public key: {} ({} ms)", pk_b, time_gen_b);

    current = time::precise_time_ns();
    let shared_a = class_group_action(&inst, pk_b, sk_a);
    let time_shared_a = (time::precise_time_ns()-current)/1000000;
    println!("Alice's shared secret: {} ({} ms)", shared_a, time_shared_a);


    current = time::precise_time_ns();
    let shared_b = class_group_action(&inst, pk_a, sk_b);
    let time_shared_b = (time::precise_time_ns()-current)/1000000;
    println!("Bob's shared secret: {} ({} ms)", shared_b, time_shared_b);

    assert_eq!(shared_a, shared_a);
    println!("Shared secret match!");
}
