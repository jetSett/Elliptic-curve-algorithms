use super::*;

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