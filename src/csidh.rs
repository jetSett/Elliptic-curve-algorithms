use rand::Rng;

use crate::finite_fields::*;
use crate::field::*;
use crate::elliptic_curves::*;
use crate::elliptic_curves::fp_elliptic_curves::*;

const N_PRIMES : usize = 13;

pub type Integer = gmp::mpz::Mpz;

type PrimeList = [Integer; N_PRIMES];

pub struct CSIDHInstance{
    pub l: PrimeList,
    pub p: Integer,
}

declare_finite_field!(K, Integer, Integer::from_str_radix("37118532150319619", 10).unwrap(), _m);

pub type PublicKey = K;

pub type SecretKey =  Vec<i32>;

pub fn check_well_defined(inst : &CSIDHInstance){
    let L = &inst.l;
    let P = &inst.p;

    let mut prod = Integer::from(1);
    for i in 0..N_PRIMES{
        prod *= L[i].clone();
    }
    assert_eq!(*P, Integer::from(4)*prod-Integer::from(1));
}

fn is_supersingular(inst : &CSIDHInstance, ell : &EllipticCurve<K>) -> bool{
    let L = &inst.l;
    let P = &inst.p;

    loop{
        let p = EllipticCurve::unsigne_point(ell.sample_point()); // We test multiple point if necessary

        let mut d : Integer = Integer::from(1);

        for i in 0..N_PRIMES{
            let li = L[i].clone();
            let qi = ell.scalar_mult_unsigned((P+1)/li.clone(), p.clone());

            if ell.scalar_mult_unsigned(li.clone(), qi.clone()) != UnsignedProjPoint::infinite_point() {
                return false;
            }
            if qi != UnsignedProjPoint::infinite_point(){
                d *= li;
            }
            if &d*&d > &Integer::from(16)*P{
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

    let p_normalized = p.clone().normalize();

    let mut order = Integer::from(2);

    // we start at i = 2 because of special doubling case
    let mut t = ell.x_dbl(p_normalized.clone()); // t = [i]p will iterate over elements of <p>
    let mut t_minus_1 = p_normalized;

    while t != UnsignedProjPoint::infinite_point(){
        order += Integer::from(1);

        let _temp = t.clone();
        t = ell.x_add(t, p.clone(), t_minus_1);
        t_minus_1 = _temp;
    }

    order
}

fn isogeny(ell : &EllipticCurve<K>, point : &UnsignedProjPoint<K>, mut q : UnsignedProjPoint<K>, k : Integer) -> Result<(EllipticCurve<K>, UnsignedProjPoint<K>), ()>{
    assert!(ell.is_montgomery());
    assert!(k>=Integer::from(3));
    assert!(&k%2 != Integer::from(0));

    let p_normalized = point.clone().normalize();
    let p = point.clone();
    q = q.normalize();


    // we start at i = 2 because of special doubling case
    let mut t = p_normalized.clone();  // t = [i]p will iterate over elements of <p>
    let mut t_minus_1 = UnsignedProjPoint::infinite_point();

    let mut pi = UnsignedProjPoint::finite_point(K::from_int(1));
    let mut sigma = K::from_int(0);
    let mut projection_x = K::from_int(1);
    let mut projection_z = K::from_int(1);

    let mut i = Integer::from(1);
    while &Integer::from(2)*&i < k{
    //     println!("{}", t);
    // while t != UnsignedProjPoint::infinite_point(){
        if t.x == K::from_int(0){ // point of order 2
            return Err(());
        }

        pi.x *= t.x.clone();
        pi.x *= t.x.clone();

        pi.z *= t.z.clone();
        pi.z *= t.z.clone();

        let s = t.x.clone()/t.z.clone() - t.z.clone()/t.x.clone();

        sigma += K::from_int(2)*s;

        let px = &t.x*&q.x - t.z.clone();
        let pz = &q.x*&t.z - t.x.clone();

        projection_x *= px.clone()*px;
        projection_z *= pz.clone()*pz;

        if i == Integer::from(1){
            let _temp = t.clone();
            t = ell.x_dbl(p.clone());
            t_minus_1 = _temp;
        }else{
            let _temp = t.clone();
            t = ell.x_add(t, p.clone(), t_minus_1);
            t_minus_1 = _temp;
        }
        i += Integer::from(1);
    }

    pi = pi.normalize();
    Ok((
        EllipticCurve::new_montgomery(pi.x*(&ell.a_2 - &(K::from_int(3)*sigma))),
        UnsignedProjPoint{
                x: q.x*projection_x,
                z: projection_z
            }
    ))

}


pub fn verify_public_key(inst : &CSIDHInstance, pk : PublicKey) -> bool{
    is_supersingular(inst, &EllipticCurve::new_montgomery(pk))
}

pub fn class_group_action(inst : &CSIDHInstance, pk : PublicKey, mut sk : SecretKey) -> PublicKey{
    let L = &inst.l;
    let P = &inst.p;

    let mut sum_abs : u32 = 0;


    for i in 0..sk.len(){
        sum_abs += (if sk[i] >= 0 { sk[i] } else {-sk[i] }) as u32;
    }

    let mut ell = EllipticCurve::new_montgomery(pk);

    while sum_abs > 0{
        let x = K::new(Integer::sample_uniform(&Integer::from(0), &(P-Integer::from(1))));
        // println!("{}", sum_abs);
        let compute_rhs = | x | {
            &(&( &(x*x)*x )*x + &((&ell.a_2)*x)*x )+ x 
        };

        let s = compute_rhs(&x).legendre_symbol() as i32;
        if s == 0{
            continue;
        }

        let uns_p = UnsignedProjPoint::finite_point(x);

        let mut s_vec = vec!();
        let mut k = Integer::from(1);

        for i in 0..sk.len(){
            if sk[i]* s > 0{
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

            let (ell_, q_point_) = match isogeny(&ell, &r_point, q_point.clone(), L[i].clone()){
                Err(()) => {
                    println!("Erreur");
                    continue;
                    },
                Ok(pair) => pair
            };

            ell = ell_;
            q_point = q_point_;

            assert!(is_supersingular(inst, &ell));
            sk[i] -= s as i32;

            k /= L[i].clone();

            sum_abs -= 1;
        }
    }
    ell.a_2
}

pub fn sample_keys(inst : &CSIDHInstance, m : i32) -> (PublicKey, SecretKey){
    let mut rng = rand::thread_rng();
    let mut sk : SecretKey = vec!();
    for _i in 0..N_PRIMES{
        sk.push(rng.gen_range(-m, m) as i32);
    }

    let pk = class_group_action(inst, PublicKey::from_int(0), sk.clone());
    (pk, sk)
}

#[cfg(test)]
mod test;