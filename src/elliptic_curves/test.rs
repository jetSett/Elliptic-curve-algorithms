use super::*;

use crate::finite_fields::*;
use crate::field::*;

const P : Integer = 10169;

declare_finite_field!(K, P, m10169);

fn sample_point(ell : &EllipticCurve<K>) -> ProjKPoint<K>{
    let mut p = ProjKPoint::FinPoint(
        K::new(0), K::new(0)
    );

    while !(ell.is_on_curve(&p)) {
        p = ProjKPoint::FinPoint(
            K::new(rand::random::<Integer>()%P),
            K::new(rand::random::<Integer>()%P),
        );
    }
    p
}

fn sample_elliptic_curve() -> EllipticCurve<K>{
    let mut ell = EllipticCurve::<K>{
                a_1: K::from_int(rand::random::<Integer>()%P),
                a_3: K::from_int(rand::random::<Integer>()%P),

                a_2: K::from_int(rand::random::<Integer>()%P),
                a_4: K::from_int(rand::random::<Integer>()%P),
                a_6: K::from_int(rand::random::<Integer>()%P),
                };
    while ell.discriminant() == K::new(0){
        ell = EllipticCurve::<K>{
                    a_1: K::from_int(rand::random::<Integer>()%P),
                    a_3: K::from_int(rand::random::<Integer>()%P),

                    a_2: K::from_int(rand::random::<Integer>()%P),
                    a_4: K::from_int(rand::random::<Integer>()%P),
                    a_6: K::from_int(rand::random::<Integer>()%P),
                    };
    }
    ell
}

#[test]
fn test_discriminant(){
    assert_eq!(EllipticCurve::<K>::new_reduced_weierstrass(K::new(0), K::new(0)).discriminant(), 
                K::new(0));

    assert_eq!(EllipticCurve::<K>{
                a_1: K::from_int(0),
                a_3: K::from_int(0),

                a_2: K::from_int(1),
                a_4: K::from_int(0),
                a_6: K::from_int(0),

                }.discriminant(), 
                K::new(0));

    //Not sure of the discr of the previous curves

    // assert_eq!(EllipticCurve::<K>::new_reduced_weierstrass(K::new(-1), K::new(0)).discriminant(), 
    //             K::new(64));

    // assert_eq!(EllipticCurve::<K>::new_reduced_weierstrass(K::new(-3), K::new(3)).discriminant(), 
    //             K::new(2160));
}

#[test]
fn reduced_weierstrass_work() {
    for _i in 1..10{
        let ell = sample_elliptic_curve();
        assert!(ell.to_reduced_weierstrass().is_reduced_weierstrass());
    }
}


#[test]
fn j_invariant_constant_weierstrass() {
    for _i in 1..10{
        let ell = sample_elliptic_curve();
        assert_eq!(ell.j_invariant(), ell.to_reduced_weierstrass().j_invariant());
    }
}

#[test]
fn point_addition_commut(){
    let ell = sample_elliptic_curve();
    for _i in 1..10{
        let p1 = sample_point(&ell);
        let p2 = sample_point(&ell);
        assert_eq!(ell.add_points(p1, p2), ell.add_points(p2, p1));
    }
}

#[test]
fn point_addition_assoc(){
    let ell = sample_elliptic_curve();
    for _i in 1..10{
        let p1 = sample_point(&ell);
        let p2 = sample_point(&ell);
        let p3 = sample_point(&ell);
        assert_eq!(ell.add_points(p1, ell.add_points(p2, p3)), ell.add_points(ell.add_points(p1, p2), p3));
    }
}

fn trivial_scalar_mult(ell : &EllipticCurve<K>, n : Integer, point : ProjKPoint<K>) -> ProjKPoint<K>{
    if n < 0 {
        return trivial_scalar_mult(ell, -n, ell.neg_point(point));
    }

    let mut p2 = ProjKPoint::InfPoint;
    for _i in 1..(n+1){
        p2 = ell.add_points(p2, point);
    }
    p2
}

#[test]
fn scalar_mult_correct() {
    let ell = sample_elliptic_curve();
    for _i in 1..10{
        let n = rand::random::<Integer>() % 100;
        let p = sample_point(&ell);
        let p_trivial = trivial_scalar_mult(&ell, n, p);
        let p_ladder = ell.scalar_mult(n, p);
        assert_eq!(p_trivial, p_ladder);
    }
}

#[test]
fn velu_isogeny_is_isogeny(){
    for _i in 1..10{
        let ell = sample_elliptic_curve().to_reduced_weierstrass();
        let p = sample_point(&ell);

        assert_eq!(ell.velu_projection(&p, ProjKPoint::InfPoint), ProjKPoint::InfPoint);
    }
}

#[test]
fn velu_isogeny_ip_kernel(){
    for _i in 1..10{
        let ell = sample_elliptic_curve().to_reduced_weierstrass();
        let p = sample_point(&ell);

        for _j in 1..10{
            let n = rand::random::<Integer>()%500;

            assert_eq!(ell.velu_projection(&p, ell.scalar_mult(n, p)), ProjKPoint::InfPoint);
        }
    }
}

#[test]
fn velu_formula_match_isogeny(){
    for _i in 1..10{
        let ell = sample_elliptic_curve().to_reduced_weierstrass();
        let p = sample_point(&ell);

        let ell_velu = ell.velu_curve(&p);

        for _j in 1..10{
            let q = sample_point(&ell);
            assert!(ell_velu.is_on_curve(&ell.velu_projection(&p, q)));
        }
    }
}