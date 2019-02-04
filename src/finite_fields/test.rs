use super::*;

declare_finite_field!(GL8001047, 8001047, m8001047);
declare_finite_field!(GL5483, 5483, m5483);
declare_finite_field!(GL1009, 1009, m1009); // = 1 mod 8
declare_finite_field!(GL2, 2, m2);

#[test]
fn addition_zero() {
    for _i in 0 .. 100{
        let a = rand::random::<Integer>();
        assert_eq!(GL2::new(a) + GL2::from_int(0), GL2::new(a));
        assert_eq!(GL5483::new(a) + GL5483::from_int(0), GL5483::new(a));
        assert_eq!(GL8001047::new(a) + GL8001047::from_int(0),GL8001047::new(a));
    }
}

#[test]
fn addition_commut() {
    for _i in 0 .. 100{
        let a = rand::random::<Integer>();
        let b = rand::random::<Integer>();
        assert_eq!(GL2::new(a) + GL2::new(b), GL2::new(b) + GL2::new(a));
        assert_eq!(GL5483::new(a) + GL5483::new(b), GL5483::new(b) + GL5483::new(a));
        assert_eq!(GL8001047::new(a) + GL8001047::new(b),GL8001047::new(b) + GL8001047::new(a));
    }
}

#[test]
fn is_square_works(){
    assert_eq!(GL8001047::new(-1).legendre_symbol(), -1);
    assert_eq!(GL5483::new(-1).legendre_symbol(), -1);
    assert_eq!(GL2::new(-1).legendre_symbol(), 1);

    assert_eq!(GL8001047::new(2).legendre_symbol(), 1);
    assert_eq!(GL5483::new(2).legendre_symbol(), -1);

    for _i in 0 .. 100{
        let a = rand::random::<Integer>();
        if a%5483 != 0{
            assert_eq!((GL5483::new(a)*GL5483::new(a)).legendre_symbol(), 1);
        }
        if a%8001047 != 0{
            assert_eq!((GL8001047::new(a)*GL8001047::new(a)).legendre_symbol(), 1);
        }

    }
}

#[test]
fn square_root_works(){

    for _i in 0 .. 100{
        let a = rand::random::<Integer>();
        let a_gl1009 = GL1009::new(a);
        let a_gl5483 = GL5483::new(a);
        let a_gl8001047 = GL8001047::new(a);

        if a_gl1009.legendre_symbol() == 1{
            let sq = a_gl1009.square_root();
            assert_eq!(sq*sq, a_gl1009);
        }

        if a_gl5483.legendre_symbol() == 1{
            let sq = a_gl5483.square_root();
            assert_eq!(sq*sq, a_gl5483);
        }

        if a_gl8001047.legendre_symbol() == 1{
            let sq = a_gl8001047.square_root();
            assert_eq!(sq*sq, a_gl8001047);
        }

        let sq8001047 = (a_gl8001047*a_gl8001047).square_root();
        let sq5483 = (a_gl5483*a_gl5483).square_root();
        let sq1009 = (a_gl1009*a_gl1009).square_root();

        assert!(sq1009 == a_gl1009 || sq1009 == -a_gl1009);
        assert!(sq5483 == a_gl5483 || sq5483 == -a_gl5483);
        assert!(sq8001047 == a_gl8001047 || sq8001047 == -a_gl8001047);
    }
}