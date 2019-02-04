use rand::Rng;

use crate::finite_fields::*;
use crate::field::*;
use crate::elliptic_curves::*;
use crate::elliptic_curves::fp_elliptic_curves::*;

const N_PRIMES : usize = 6;

const L : [Integer; N_PRIMES] = [3, 5, 7, 11, 13, 17];

const P : Integer = 1021019; // NB : P = 4*L_0*...*L_{N_PRIMES-1} -1

declare_finite_field!(K, P, _m);

pub type PublicKey = K;

pub fn check_well_defined(){
    let mut prod = 1;
    for i in 0..N_PRIMES{
        prod *= L[i];
    }
    assert_eq!(P, 4*prod-1);
}

fn is_supersingular(ell : &EllipticCurve<K>) -> bool{
    loop{
        let p = EllipticCurve::unsigne_point(ell.sample_point()); // We test multiple point if necessary

        let mut d : Integer = 1;

        for i in 0..N_PRIMES{
            let li = L[i];
            let qi = ell.scalar_mult_unsigned((P+1)/li, p);

            if ell.scalar_mult_unsigned(li, qi) != UnsignedProjPoint::infinite_point() {
                return false;
            }
            if qi != UnsignedProjPoint::infinite_point(){
                d *= li;
            }
            if d*d > 16*P{
                return true;
            }
        }
    }
}

fn order_naive(ell : &EllipticCurve<K>, p : &UnsignedProjPoint<K>) -> u32{
    assert!(ell.is_montgomery());

    let p_normalized = (*p).normalize();


    let mut order : u32 = 2;

    // we start at i = 2 because of special doubling case
    let mut t = ell.x_dbl(p_normalized); // t = [i]p will iterate over elements of <p>
    // println!("{}", t);

    let mut t_minus_1 = p_normalized;

    while t != UnsignedProjPoint::infinite_point(){
        order += 1;

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

    let mut projection_x = (q.x*p_normalized.x - K::from_int(1));
    let mut projection_z = (q.x - p_normalized.x);

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

    let mut t_minus_1 = p_normalized;

    let mut pi = p_normalized.x;
    let mut sigma = p_normalized.x - K::from_int(1)/p_normalized.x;

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

pub fn class_group_action(pk : PublicKey, mut e : Vec<i8>) -> PublicKey{
    let mut rng = rand::thread_rng();
    let mut sum_abs : u32 = 0;

    let mut a = pk;

    for i in 0..e.len(){
        sum_abs += (if e[i] >= 0 { e[i] } else {-e[i] }) as u32;
    }

    let mut ell = EllipticCurve::new_montgomery(pk);

    while sum_abs > 0{
        let x = K::new(rng.gen_range(0, P-1));
        let s = (x*x*x + pk*x*x + x).legendre_symbol();

        let uns_p = UnsignedProjPoint::finite_point(x);

        let mut s_vec = vec!();
        let mut k = 1;

        for i in 0..e.len(){
            if e[i]* s > 0{
                s_vec.push(i);
                k *= L[i];
            }
        }

        if s_vec.is_empty(){
            continue;
        }

        let mut q_point = ell.scalar_mult_unsigned((P+1)/k, uns_p);

        for j in 0..s_vec.len(){
            let i = s_vec[j];

            let r_point = ell.scalar_mult_unsigned(k/L[i], q_point);

            if r_point == UnsignedProjPoint::infinite_point(){
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

            k /= L[i];
            e[i] -= s;
            sum_abs -= 1;

            a = ell.a_2;
        }
    }
    a
}

#[cfg(test)]
mod test;