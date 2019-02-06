use super::*;

fn sample_montgomery() -> EllipticCurve<K>{

    let inst = CSIDHInstance{ 
        p: Integer::from_str_radix("37118532150319619", 10).unwrap(),
        l:[Integer::from(3),Integer::from(5),Integer::from(7),
                Integer::from(11),Integer::from(13),Integer::from(17),
                Integer::from(19),Integer::from(23),Integer::from(29),
                Integer::from(31),Integer::from(37),Integer::from(41),
                Integer::from(61)]
    };

    let mut a = K::new(Integer::sample_uniform(&Integer::from(0), &Integer::from(inst.p.clone())));
    let mut ell = EllipticCurve::new_montgomery(a.clone());
    while ! (is_supersingular(&inst, &ell) && 
            ell.discriminant() != K::from_int(0)) {
        a = K::new(Integer::sample_uniform(&Integer::from(0), &Integer::from(inst.p.clone())));
        ell = EllipticCurve::new_montgomery(a);
    }
    ell
}

#[test]
fn everything_ok() {
    let inst = CSIDHInstance{ 
        p: Integer::from_str_radix("37118532150319619", 10).unwrap(),
        l:[Integer::from(3),Integer::from(5),Integer::from(7),
                Integer::from(11),Integer::from(13),Integer::from(17),
                Integer::from(19),Integer::from(23),Integer::from(29),
                Integer::from(31),Integer::from(37),Integer::from(41),
                Integer::from(61)]
    };

    check_well_defined(&inst);
}

#[test]
fn known_supersingular_ec(){

    let inst = CSIDHInstance{ 
        p: Integer::from_str_radix("37118532150319619", 10).unwrap(),
        l:[Integer::from(3),Integer::from(5),Integer::from(7),
                Integer::from(11),Integer::from(13),Integer::from(17),
                Integer::from(19),Integer::from(23),Integer::from(29),
                Integer::from(31),Integer::from(37),Integer::from(41),
                Integer::from(61)]
    };

    // Note P = 3 mod 4, 2 mod 3
    assert!(is_supersingular(&inst, &EllipticCurve::new_reduced_weierstrass(K::from_int(1), K::from_int(0)))); // supersingular iff P = 3 mod 4
    assert!(is_supersingular(&inst, &EllipticCurve::new_reduced_weierstrass(K::from_int(0), K::from_int(1)))); // supersingular iff P = 2 mod 3
}

#[test]
fn class_group_action_keep_valid() {
    let inst = CSIDHInstance{ 
        p: Integer::from_str_radix("37118532150319619", 10).unwrap(),
        l:[Integer::from(3),Integer::from(5),Integer::from(7),
                Integer::from(11),Integer::from(13),Integer::from(17),
                Integer::from(19),Integer::from(23),Integer::from(29),
                Integer::from(31),Integer::from(37),Integer::from(41),
                Integer::from(61)]
    };

    for _i in 0..10{
        let mut pk = PublicKey::from_int(0);
    
        let mut action : SecretKey = vec![0, 0, 0, 0, 0, 0];
        for j in 0..6{
            action[j] = rand::random::<i32>()%3;
        }
        assert!(verify_public_key(&inst, class_group_action(&inst, pk, action)));
    }

}

#[test]
fn class_group_action_commute(){
    let inst = CSIDHInstance{ 
        p: Integer::from_str_radix("37118532150319619", 10).unwrap(),
        l:[Integer::from(3),Integer::from(5),Integer::from(7),
                Integer::from(11),Integer::from(13),Integer::from(17),
                Integer::from(19),Integer::from(23),Integer::from(29),
                Integer::from(31),Integer::from(37),Integer::from(41),
                Integer::from(61)]
    };

    for _i in 0..10{
        let mut pk = K::from_int(0);
    
        let mut action1 : SecretKey = vec!();
        let mut action2 : SecretKey = vec!();

        for j in 0..N_PRIMES{
            action1.push(rand::random::<i32>()%2);
            action2.push(rand::random::<i32>()%2);
        }


        let pk1 = class_group_action(&inst, pk.clone(), action1.clone());
        let pk12 = class_group_action(&inst, pk1.clone(), action2.clone());


        let pk2 = class_group_action(&inst, pk, action2.clone());
        let pk21 = class_group_action(&inst, pk2, action1.clone());

        assert_eq!(pk12, pk21);
    }
}