use super::*;

use crate::finite_fields::*;

declare_finite_field!(K, 10169, m10169);

fn sample_point(ell : &EllipticCurve<K>) -> ProjKPoint<K>{
    let mut p = ProjKPoint::FinPoint(KPoint::<K>{
        x: K::new(0), y: K::new(0)
    });

    while !(ell.is_on_curve(&p)) {
        p = ProjKPoint::FinPoint(KPoint{
            x: K::new(rand::random::<Integer>()),
            y: K::new(rand::random::<Integer>()),
        });
    }
    p
}

#[test]
fn point_addition_commut(){
    let ell = EllipticCurve::<K>::new(K::new(1), K::new(5841));
    for _i in 1..10{
        let n = rand::random::<i32>() % 100;
        let p1 = sample_point(&ell);
        let p2 = sample_point(&ell);
        assert_eq!(ell.add_points(p1, p2), ell.add_points(p2, p1));
    }
}

#[test]
fn point_addition_assoc(){
    let ell = EllipticCurve::<K>::new(K::new(1), K::new(5841));
    for _i in 1..10{
        let n = rand::random::<i32>() % 100;
        let p1 = sample_point(&ell);
        let p2 = sample_point(&ell);
        let p3 = sample_point(&ell);
        assert_eq!(ell.add_points(p1, ell.add_points(p2, p3)), ell.add_points(ell.add_points(p1, p2), p3));
    }
}

fn trivial_scalar_mult(ell : &EllipticCurve<K>, point : ProjKPoint<K>, n : i32) -> ProjKPoint<K>{
    let mut p = point;
    let mut m = n;
    if n < 0 {
        m = -m;
        p = ell.neg_point(p);
    }
    let mut p2 = ProjKPoint::InfPoint;
    for _i in 1..m{
        p2 = ell.add_points(p2, p);
    }
    p2
}

#[test]
fn scalar_mult_correct() {
    let ell = EllipticCurve::<K>::new(K::new(1), K::new(5841));
    for _i in 1..10{
        let n = rand::random::<i32>() % 100;
        let p = sample_point(&ell);
        assert_eq!(ell.scalar_mult(n, p), trivial_scalar_mult(&ell, p, n));
    }
}
