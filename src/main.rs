#![feature(type_alias_enum_variants)]

use time;

#[macro_use] pub mod finite_fields;
pub mod elliptic_curves;
pub mod field;

mod csidh;

use csidh::*;

fn main() {
    let gen_param = 2;

    let mut current : u64;

    let inst = CSIDHInstance{ 
        p: Integer::from_str_radix("37118532150319619", 10).unwrap(),
        l:[Integer::from(3),Integer::from(5),Integer::from(7),
                Integer::from(11),Integer::from(13),Integer::from(17),
                Integer::from(19),Integer::from(23),Integer::from(29),
                Integer::from(31),Integer::from(37),Integer::from(41),
                Integer::from(61)]
    };

    check_well_defined(&inst);

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

    assert_eq!(shared_a, shared_b);
    println!("Shared secret match!");
}
