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
fn addition_coincide() {
    for _j in 1..10{
        let ell = sample_montgomery();
        for _i in 1..50{
            let p1 = ell.sample_point();
            let p2 = ell.sample_point();

            let p1compr = ell.compress_point(p1);
            let p2compr = ell.compress_point(p2);

            assert_eq!(ell.add_points(p1, p2), 
                    ell.decompress_point(ell.add_compressed_point(p1compr, p2compr)));

        }
    }
}

#[test]
fn scalar_multiplication_coincide(){
    for _j in 1..10{
        let ell = sample_montgomery();
        for _i in 1..50{
            let n = rand::random::<Integer>()%1000;
            let p = ell.sample_point();

            let pcompr = ell.compress_point(p);

            assert_eq!(ell.scalar_mult(n, p), 
                    ell.decompress_point(ell.scalar_mult_compressed(n, pcompr)));

        }
    }

}