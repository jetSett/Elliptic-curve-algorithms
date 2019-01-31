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

            println!("{}", d);
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

#[cfg(test)]
mod test;