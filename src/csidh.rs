use rand::Rng;

use crate::finite_fields::*;
use crate::field::*;
use crate::elliptic_curves::*;
use crate::elliptic_curves::fp_elliptic_curves::*;

const N_PRIMES : usize = 6;

type Integer = gmp::mpz::Mpz;

const L : [Integer; N_PRIMES] = [Integer::from(3), Integer::from(5), Integer::from(7), Integer::from(11), Integer::from(13), Integer::from(17)];

const P : Integer = Integer::from(1021019); // NB : P = 4*L_0*...*L_{N_PRIMES-1} -1

declare_finite_field!(K, Integer, P, _m);

pub type PublicKey = K;

pub type SecretKey =  [i32; N_PRIMES];

pub fn check_well_defined(){
    let mut prod = Integer::from(1);
    for i in 0..N_PRIMES{
        prod *= L[i];
    }
    assert_eq!(P, Integer::from(4)*prod-Integer::from(1));
}

fn is_supersingular(ell : &EllipticCurve<K>) -> bool{
    loop{
        let p = EllipticCurve::unsigne_point(ell.sample_point()); // We test multiple point if necessary

        let mut d : Integer = Integer::from(1);

        for i in 0..N_PRIMES{
            let li = L[i];
            let qi = ell.scalar_mult_unsigned((P+1)/li, p);

            if ell.scalar_mult_unsigned(li, qi) != UnsignedProjPoint::infinite_point() {
                return false;
            }
            if qi != UnsignedProjPoint::infinite_point(){
                d *= li;
            }
            if d*d > Integer::from(16)*P{
                return true;
            }
        }
    }
}

fn order_naive(ell : &EllipticCurve<K>, p : &UnsignedProjPoint<K>) -> Integer{
    assert!(ell.is_montgomery());

    if p == &UnsignedProjPoint::infinite_point(){
        return Integer::from(1);
    }

    let p_normalized = (*p).normalize();


    let mut order = Integer::from(2);

    // we start at i = 2 because of special doubling case
    let mut t = ell.x_dbl(p_normalized); // t = [i]p will iterate over elements of <p>
    // println!("{}", t);

    let mut t_minus_1 = p_normalized;

    while t != UnsignedProjPoint::infinite_point(){
        order += Integer::from(1);

        let _temp = t;
        t = ell.x_add(t, *p, t_minus_1);
        t_minus_1 = _temp;
    }

    order
}


fn velu_projection_montgomery(ell : &EllipticCurve<K>, p : &UnsignedProjPoint<K>, mut q : UnsignedProjPoint<K>) -> UnsignedProjPoint<K>{
    assert!(ell.is_montgomery());

    let p_normalized = (*p).normalize();
    q = q.normalize();

    let mut projection_x = q.x*p_normalized.x - K::from_int(1);
    let mut projection_z = q.x - p_normalized.x;

    // we start at i = 2 because of special doubling case
    let mut t = ell.x_dbl(p_normalized); // t = [i]p will iterate over elements of <p>

    let mut t_minus_1 = p_normalized;


    while t != UnsignedProjPoint::infinite_point(){
        t = t.normalize();
        projection_x *= t.x*q.x - K::from_int(1);
        projection_z *= q.x - t.x;

        let _temp = t;
        t = ell.x_add(t, *p, t_minus_1);
        t_minus_1 = _temp;
    }

    UnsignedProjPoint{
        x: q.x*projection_x,
        z: projection_z
    }
}

fn velu_formula_montgomery(ell : &EllipticCurve<K>, point : &UnsignedProjPoint<K>) -> Result<EllipticCurve<K>, ()>{
    assert!(ell.is_montgomery());
    let a = ell.a_2;

    if point.x == K::from_int(0){
        return Err(());
    }

    let p_normalized = (*point).normalize();

    // we start at i = 2 because of special doubling case
    let mut t = ell.x_dbl(p_normalized); // t = [i]p will iterate over elements of <p>

    let mut t_minus_1 = p_normalized.clone();

    let mut pi = p_normalized.x.clone();
    let mut sigma = p_normalized.x.clone() - K::from_int(1)/p_normalized.x.clone();

    while t != UnsignedProjPoint::infinite_point(){
        if t.x == K::from_int(0){
            return Err(());
        }

        t = t.normalize();
        t_minus_1 =t_minus_1.normalize();
        
        pi *= t.x;
        sigma += t.x - K::from_int(1)/t.x;

        let _temp = t;
        t = ell.x_add(t, p_normalized, t_minus_1);
        t_minus_1 = _temp;
    }

    Ok(EllipticCurve::new_montgomery(pi*(a - K::from_int(3)*sigma)))

}

pub fn verify_public_key(pk : PublicKey) -> bool{
    is_supersingular(&EllipticCurve::new_montgomery(pk))
}

pub fn class_group_action(pk : PublicKey, mut sk : SecretKey) -> PublicKey{
    let mut sum_abs : u32 = 0;


    for i in 0..sk.len(){
        sum_abs += (if sk[i] >= 0 { sk[i] } else {-sk[i] }) as u32;
    }

    let mut ell = EllipticCurve::new_montgomery(pk);

    while sum_abs > 0{
        let x = K::new(Integer::sample_uniform(&Integer::from(0), &(P-Integer::from(1))));

        let compute_rhs = | x | {
            &(&( &(x*x)*x )*x + &((&ell.a_2)*x)*x )+ x 
        };

        let s = compute_rhs(&x).legendre_symbol();

        if s == 0{
            continue;
        }

        let uns_p = UnsignedProjPoint::finite_point(x);

        let mut s_vec = vec!();
        let mut k = Integer::from(1);

        for i in 0..sk.len(){
            if sk[i]* (s as i32) > 0{
                s_vec.push(i);
                k *= L[i].clone();
            }
        }

        if s_vec.is_empty(){
            continue;
        }

        let mut q_point = ell.scalar_mult_unsigned((P.clone()+Integer::from(1))/k.clone(), uns_p);

        for j in 0..s_vec.len(){
            let i = s_vec[s_vec.len()-j-1];

            let r_point = ell.scalar_mult_unsigned(k.clone()/L[i].clone(), q_point.clone());

            if &r_point == &UnsignedProjPoint::infinite_point(){
                continue;
            }

            q_point = velu_projection_montgomery(&ell, &r_point, q_point);
            ell = match velu_formula_montgomery(&ell, &r_point){
                Err(()) => {
                    continue;
                    },
                Ok(ell_) => ell_,
            };

            assert!(is_supersingular(&ell));
            sk[i] -= s as i32;

            k /= L[i].clone();

            sum_abs -= 1;
        }
    }
    ell.a_2
}

pub fn sample_keys(m : i32) -> (PublicKey, SecretKey){
    let mut rng = rand::thread_rng();
    let mut sk : SecretKey = [0, 0, 0, 0, 0, 0];
    for i in 0..N_PRIMES{
        sk[i] = rng.gen_range(-m, m) as i32;
    }
    let pk = class_group_action(PublicKey::from_int(0), sk);
    (pk, sk)
}

#[cfg(test)]
mod test;