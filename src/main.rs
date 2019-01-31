#![feature(type_alias_enum_variants)]

#[macro_use] pub mod finite_fields;
pub mod elliptic_curves;
pub mod field;

use finite_fields::*;
use elliptic_curves::*;

use field::*;

const P : Integer = 10169;

declare_finite_field!(K, P, m10169);

fn sample_point(ell : &EllipticCurve<K>) -> ProjKPoint<K>{
    let mut p = ProjKPoint::FinPoint(
        K::new(0), K::new(0)
    );

    while !(ell.is_on_curve(&p)) {
        p = ProjKPoint::FinPoint(
            K::new(rand::random::<Integer>()),
            K::new(rand::random::<Integer>()),
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

fn trivial_order(ell: &EllipticCurve<K>, p : &ProjKPoint<K>) -> u32{
    let mut p1 = *p;

    let mut i = 1;

    while p1 != ProjKPoint::InfPoint && ell.is_on_curve(&p1){
        // println!("{}", p1);
        p1 = ell.add_points(p1, *p);
        i+=1;
    }

    i
}

fn main() {

    let ell = sample_elliptic_curve().to_reduced_weierstrass();

    let p = sample_point(&ell);

    println!("Found ! computing order...\n");

    let orderp = trivial_order(&ell, &p);

    println!("\nOrder of {} : {}", p, orderp);


    let q = sample_point(&ell);
    let orderq = trivial_order(&ell, &q);

    println!("\nOrder of {} : {}", q, orderq);


    println!("{} --> {} ({})", q, ell.velu_projection(&p, q), orderp%orderq);
}
