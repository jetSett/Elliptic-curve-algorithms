use super::*;

declare_finite_field!(GL8001047, 8001047, m1);
declare_finite_field!(GL5483, 5483, m2);
declare_finite_field!(GL2, 2, m3);

#[test]
fn addition_zero() {
    for _i in 0 .. 100{
        let a = rand::random::<Integer>();
        assert_eq!(GL2::new(a) + GL2::zero(), GL2::new(a));
        assert_eq!(GL5483::new(a) + GL5483::zero(), GL5483::new(a));
        assert_eq!(GL8001047::new(a) + GL8001047::zero(),GL8001047::new(a));
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