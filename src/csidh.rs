use rand::Rng;

use crate::finite_fields::*;
use crate::field::*;
use crate::elliptic_curves::*;

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
        let p = ell.efficient_sample_point(); // We test multiple point if necessary

        let mut d : Integer = 1;

        for i in 0..N_PRIMES{
            let li = L[i];
            let qi = ell.scalar_mult((P+1)/li, p);

            if ell.scalar_mult(li, qi) != ProjKPoint::InfPoint {
                return false;
            }
            if qi != ProjKPoint::InfPoint{
                d *= li;
            }
            if d*d > 16*P{
                return true;
            }
        }
    }
}


pub fn verify_public_key(pk : PublicKey) -> bool{
    is_supersingular(&EllipticCurve::new_montgomery(pk))
}

pub fn class_group_action(pk : K, mut e : Vec<i8>){
    let mut rng = rand::thread_rng();
    let mut sum_abs = 0;
    for i in 0..e.len(){
        sum_abs += if e[i] >= 0 { e[i] } else {-e[i] };
    }

    let ell = EllipticCurve::new_montgomery(pk);
    println!("{}", ell);

    let random_point = ell.efficient_sample_point();

    while sum_abs > 0{
        let x = K::new(rng.gen_range(0, P-1));
        let s = (x*x*x + pk*x*x + x).legendre_symbol();
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

        let mut q_point = ell.scalar_mult((P+1)/k, random_point);
        for j in 0..s_vec.len(){
            let i = s_vec[j];

            let r_point = ell.scalar_mult(k/L[i], q_point);
            let ell2 = ell.velu_curve(&r_point);
            assert!(is_supersingular(&ell2));

            
            let q_point = ell.velu_projection(&r_point, q_point);

            k /= L[i];
            e[i] -= s;
            sum_abs -= 1;
            println!("{}", ell2);
        }
        break;
    }
}

#[cfg(test)]
mod test;