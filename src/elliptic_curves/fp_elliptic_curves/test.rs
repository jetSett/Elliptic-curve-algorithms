use super::*;

const P : Integer = 10169;

declare_finite_field!(K, P, m10169);

fn sample_montgomery() -> EllipticCurve<K>{
    let mut ell = EllipticCurve::new_montgomery(K::new(rand::random::<Integer>()%P));
    while ell.discriminant() == K::new(0){
        ell = EllipticCurve::new_montgomery(K::new(rand::random::<Integer>()%P));
    }
    ell
}

#[test]
fn x_double_coincide(){
    for _j in 1..10{
        let ell = sample_montgomery();
        for _i in 1..50{
            let p = ell.sample_point();
            let p2 = ell.add_points(p, p);
            let p_int = EllipticCurve::unsigne_point(p);

            assert_eq!(ell.x_dbl(p_int), 
                       EllipticCurve::unsigne_point(p2));
        }
    }
}

#[test]
fn x_add_coincide(){
    for _j in 1..10{
        let ell = sample_montgomery();
        for _i in 1..50{
            let p1 = ell.sample_point();
            let p2 = ell.sample_point();
            let p_plus = ell.add_points(p1, p2);
            let p_moins = ell.add_points(p1, ell.neg_point(p2));

            let p1_uns = EllipticCurve::unsigne_point(p1); 
            let p2_uns = EllipticCurve::unsigne_point(p2); 
            let p_plus_uns = EllipticCurve::unsigne_point(p_plus); 
            let p_moins_uns = EllipticCurve::unsigne_point(p_moins); 

            assert_eq!(ell.x_add(p1_uns, p2_uns, p_moins_uns), 
                       p_plus_uns);
        }
    }
}

#[test]
fn scalar_multiplication_coincide_scalar(){
    for _j in 1..10{
        let ell = sample_montgomery();
        for _i in 1..50{
            let n = rand::random::<Integer>()%1000;
            let p = ell.sample_point();

            let np = ell.scalar_mult(n, p);

            let unsign_p = EllipticCurve::unsigne_point(p);

            assert_eq!(EllipticCurve::unsigne_point(np), 
                    ell.scalar_mult_unsigned(n, unsign_p));

        }
    }
}