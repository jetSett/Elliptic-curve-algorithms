use super::*;

fn sample_montgomery() -> EllipticCurve<K>{
    let mut a = K::from_int(rand::random::<Integer>());
    while ! (is_supersingular(&EllipticCurve::new_montgomery(a)) && EllipticCurve::new_montgomery(a).discriminant() != K::from_int(0)) {
        a = K::from_int(rand::random::<Integer>());
    }
    EllipticCurve::new_montgomery(a)
}

#[test]
fn everything_ok() {
    check_well_defined();
}

#[test]
fn known_supersingular_ec(){
    // Note P = 3 mod 4, 2 mod 3
    assert!(is_supersingular(&EllipticCurve::new_reduced_weierstrass(K::new(1), K::new(0)))); // supersingular iff P = 3 mod 4
    assert!(is_supersingular(&EllipticCurve::new_reduced_weierstrass(K::new(0), K::new(1)))); // supersingular iff P = 2 mod 3
}

#[test]
fn class_group_action_keep_valid() {
    for _i in 0..10{
        let mut pk = PublicKey::new(0);
    
        let mut action = vec![0, 0, 0, 0, 0, 0];
        for j in 0..6{
            action[j] = rand::random::<i8>()%3;
        }
        assert!(verify_public_key(class_group_action(pk, action)));
    }

}