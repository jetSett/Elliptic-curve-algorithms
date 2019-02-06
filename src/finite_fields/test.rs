use super::*;

pub type Integer = gmp::mpz::Mpz;

declare_finite_field!(GL8001047, Integer, Integer::from(8001047), m8001047);
declare_finite_field!(GL5483, Integer, Integer::from(5483), m5483);
declare_finite_field!(GL1009, Integer, Integer::from(1009), m1009); // = 1 mod 8
declare_finite_field!(GL2, Integer, Integer::from(2), m2);

#[test]
fn addition_zero() {
    for _i in 0 .. 100{
        let a = Integer::sample_uniform(&Integer::from(0), &Integer::from(8001047));
        assert_eq!(GL2::new(a.clone()) + GL2::from_int(0), GL2::new(a.clone()));
        assert_eq!(GL5483::new(a.clone()) + GL5483::from_int(0), GL5483::new(a.clone()));
        assert_eq!(GL8001047::new(a.clone()) + GL8001047::from_int(0),GL8001047::new(a.clone()));
    }
}

#[test]
fn addition_commut() {
    for _i in 0 .. 100{
        let a = Integer::sample_uniform(&Integer::from(0), &Integer::from(8001047));
        let b = Integer::sample_uniform(&Integer::from(0), &Integer::from(8001047));
        assert_eq!(GL2::new(a.clone()) + GL2::new(b.clone()), GL2::new(b.clone()) + GL2::new(a.clone()));
        assert_eq!(GL5483::new(a.clone()) + GL5483::new(b.clone()), GL5483::new(b.clone()) + GL5483::new(a.clone()));
        assert_eq!(GL8001047::new(a.clone()) + GL8001047::new(b.clone()),GL8001047::new(b.clone()) + GL8001047::new(a.clone()));
    }
}

#[test]
fn is_square_works(){
    assert_eq!(GL8001047::new(Integer::from(-1)).legendre_symbol(), -1);
    assert_eq!(GL5483::new(Integer::from(-1)).legendre_symbol(), -1);
    assert_eq!(GL2::new(Integer::from(-1)).legendre_symbol(), 1);

    assert_eq!(GL8001047::new(Integer::from(2)).legendre_symbol(), 1);
    assert_eq!(GL5483::new(Integer::from(2)).legendre_symbol(), -1);

    for _i in 0 .. 100{
        let a = Integer::sample_uniform(&Integer::from(0), &Integer::from(8001047));
        if a.clone()%5483 != Integer::from(0){
            assert_eq!((GL5483::new(a.clone())*GL5483::new(a.clone())).legendre_symbol(), 1);
        }
        if a.clone()%8001047 != Integer::from(0){
            assert_eq!((GL8001047::new(a.clone())*GL8001047::new(a.clone())).legendre_symbol(), 1);
        }

    }
}

#[test]
fn square_root_works(){

    for _i in 0 .. 100{
        let a = Integer::sample_uniform(&Integer::from(0), &Integer::from(8001047));
        let a_gl1009 = GL1009::new(a.clone());
        let a_gl5483 = GL5483::new(a.clone());
        let a_gl8001047 = GL8001047::new(a);

        if a_gl1009.legendre_symbol() == 1{
            let sq = a_gl1009.clone().square_root();
            assert_eq!(sq.clone()*sq, a_gl1009);
        }

        if a_gl5483.legendre_symbol() == 1{
            let sq = a_gl5483.clone().square_root();
            assert_eq!(sq.clone()*sq, a_gl5483);
        }

        if a_gl8001047.legendre_symbol() == 1{
            let sq = a_gl8001047.clone().square_root();
            assert_eq!(sq.clone()*sq, a_gl8001047);
        }

        let sq8001047 = (a_gl8001047.clone()*a_gl8001047.clone()).square_root();
        let sq5483 = (a_gl5483.clone()*a_gl5483.clone()).square_root();
        let sq1009 = (a_gl1009.clone()*a_gl1009.clone()).square_root();

        assert!(sq1009 == a_gl1009 || sq1009 == -a_gl1009);
        assert!(sq5483 == a_gl5483 || sq5483 == -a_gl5483);
        assert!(sq8001047 == a_gl8001047 || sq8001047 == -a_gl8001047);
    }
}