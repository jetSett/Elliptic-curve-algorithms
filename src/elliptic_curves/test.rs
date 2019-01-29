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
        let p1 = sample_point(&ell);
        let p2 = sample_point(&ell);
        assert_eq!(ell.add_points(p1, p2), ell.add_points(p2, p1));
    }
}

#[test]
fn point_addition_assoc(){
    let ell = EllipticCurve::<K>::new(K::new(1), K::new(5841));
    for _i in 1..50{
        let p1 = sample_point(&ell);
        let p2 = sample_point(&ell);
        let p3 = sample_point(&ell);
        assert_eq!(ell.add_points(p1, ell.add_points(p2, p3)), ell.add_points(ell.add_points(p1, p2), p3));
    }
}

fn trivial_scalar_mult(ell : &EllipticCurve<K>, n : i32, point : ProjKPoint<K>) -> ProjKPoint<K>{
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
    let ell = EllipticCurve::<K>::new(K::new(1), K::new(5841));
    for _i in 1..50{
        let n = rand::random::<i32>() % 100;
        let p = sample_point(&ell);
        let p_trivial = trivial_scalar_mult(&ell, n, p);
        let p_ladder = ell.scalar_mult(n, p);
        assert_eq!(p_trivial, p_ladder);
    }
}
