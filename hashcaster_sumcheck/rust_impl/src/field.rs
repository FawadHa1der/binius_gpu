use std::{alloc::System, ops::{Add, AddAssign, BitAnd, BitAndAssign, Mul, MulAssign}};
use crate::{backend::autodetect::mul_128, precompute::{cobasis_frobenius_table::COBASIS_FROBENIUS, cobasis_table::COBASIS, frobenius_table::FROBENIUS}, utils::{u128_rand, binary128_to_bits}};
use bytemuck::{AnyBitPattern, NoUninit, Pod, Zeroable};
use num_traits::{One, Zero};

use rand::Rng;

use binius_field::{arithmetic_traits::TaggedSquare, BinaryField128b, BinaryField128bPolyval, ExtensionField};
use binius_field::{
	arithmetic_traits::{Broadcast, InvertOrZero, MulAlpha, Square},
	underlier::{NumCast, UnderlierType, UnderlierWithBitOps, WithUnderlier, U1, U2, U4},
	BinaryField, PackedField, Field, BinaryField1b
};
use binius_math::Matrix;


#[repr(transparent)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Pod, Zeroable)]
pub struct F128 {
    // pub(crate) raw: u128,
    pub inner_binius_field: BinaryField128b
}

// pub trait F128{
//     fn new(x: bool) -> Self;
//     fn from_raw(raw: u128) -> Self;
//     fn raw(&self) -> u128;
//     fn into_raw(self) -> u128;
//     fn rand<RNG: Rng>(rng: &mut RNG) -> Self;
//     fn frob(&self, k: i32) -> Self;
//     fn basis(i: usize) -> Self;
//     fn cobasis(i: usize) -> Self;
// }

impl F128 {
    pub fn new(x: bool) -> Self {
        if x {Self::one()} else {Self::zero()}
    }

    pub fn from_raw(raw: u128) -> Self {
        Self {
            inner_binius_field: BinaryField128b::new(raw),
        }
    }

    pub fn from_binius_field(binius_field: BinaryField128b) -> Self {
        Self {
            inner_binius_field: binius_field,
        }
    }

    pub fn raw(&self) -> u128 {
        u128::from(self.inner_binius_field.val())
    }

    pub fn into_raw(self) -> u128 {
        self.inner_binius_field.val()
    }

    pub fn rand<RNG: Rng>(rng: &mut RNG) -> Self {
        Self::from_raw(u128_rand(rng))
    }
    pub fn frob_direct(&self, k: u32) -> Self {
        let mut result = *self;
        for _ in 0..k {
            result = result.square(); // Implement square() to perform \(\alpha \mapsto \alpha^2\)
        }
        result
    }

    /// This function is not efficient.
    pub fn frob(&self, mut k: i32) -> Self {
        if k < 0 {
            k *= -1;
            k %= 128;
            k *= -1;
            k += 128;
        } else {
            k %= 128
        }
        let matrix = &FROBENIUS[k as usize]; 
        let mut ret = 0;
        let vec_bits = binary128_to_bits(self.inner_binius_field);
        for i in 0..128 {
            if vec_bits[i] {ret ^= matrix[i]}
        }
        F128::from_raw(ret)
    }

    pub fn basis(i: usize) -> Self {
        assert!(i < 128);
        let f =  <BinaryField128b as ExtensionField<BinaryField1b>>::basis(i).unwrap();
        // let polyval = BinaryField128bPolyval::from(f);
        
        // let polyval = <BinaryField128bPolyval as ExtensionField<BinaryField1b>>::basis(i);
        // let iterator = <BinaryField128b as ExtensionField<BinaryField1b>>::iter_bases(&f);
        // for elem in iterator {
        //     println!("{:?}", elem);
        // }
        
        // let transformed_polyval = BinaryField128b::from(polyval.unwrap());
        // let unwraped_f = f.unwrap();
        // assert!(transformed_polyval == unwraped_f);
        
        // let f: Result<BinaryField128b, binius_field::Error> =  BinaryField128b::basis(i);
        // let f: BinaryField128b = f.unwrap();
        // println!("f: {:?}", f);

        Self::from_binius_field(f)
        // Self::from_raw(1 << i)

    }

    pub fn cobasis(i: usize) -> Self {
        Self::from_raw(COBASIS[i])
    }
}

impl Zero for F128 {
    fn zero() -> Self {
        Self{inner_binius_field:BinaryField128b::ZERO}
    }

    fn is_zero(&self) -> bool {
        self.inner_binius_field == BinaryField128b::ZERO
    }
}

impl One for F128 {
    fn one() -> Self {
        Self{inner_binius_field: BinaryField128b::ONE}
    }
}

impl Add<F128> for F128 {
    type Output = F128;

    fn add(self, rhs: Self) -> Self::Output {
        Self{inner_binius_field: self.inner_binius_field + rhs.inner_binius_field}
    }
}

impl Add<&F128> for F128 {
    type Output = F128;

    fn add(self, rhs: &F128) -> Self::Output {
        Self{inner_binius_field: self.inner_binius_field + rhs.inner_binius_field}
    }
}

impl BitAnd<F128> for F128 {
    type Output = F128;
    
    fn bitand(self, rhs: F128) -> Self::Output {
        Self{inner_binius_field: BinaryField128b::new( self.into_raw() & rhs.into_raw() )}
    }
}

impl BitAnd<&F128> for F128 {
    type Output = F128;
    
    fn bitand(self, rhs: &F128) -> Self::Output {
        Self{inner_binius_field: BinaryField128b::new(self.into_raw() & rhs.raw())}
    }
}

impl AddAssign<F128> for F128 {
    fn add_assign(&mut self, rhs: F128) {
        self.inner_binius_field += rhs.inner_binius_field
    }
}

impl AddAssign<&F128> for F128 {
    fn add_assign(&mut self, rhs: &F128) {
        self.inner_binius_field += rhs.inner_binius_field
    }
}

// impl BitAndAssign<F128> for F128 {
//     fn bitand_assign(&mut self, rhs: F128) {
//         self.inner_binius_field &= rhs.inner_binius_field;
//     }
// }

// impl BitAndAssign<&F128> for F128 {
//     fn bitand_assign(&mut self, rhs: &F128) {
//         self.raw &= rhs.raw();
//     }
// }

impl Mul<F128> for F128 {
    type Output = F128;

    fn mul(self, rhs: F128) -> Self::Output {
        Self{inner_binius_field: self.inner_binius_field * rhs.inner_binius_field}
    }
}

impl Mul<&F128> for F128 {
    type Output = F128;

    fn mul(self, rhs: &F128) -> Self::Output {
        // Self::from_raw(mul_128(self.into_raw(), rhs.raw()))
        Self{inner_binius_field: self.inner_binius_field * rhs.inner_binius_field}
    }
}

impl MulAssign<F128> for F128 {
    fn mul_assign(&mut self, rhs: F128) {
        *self = *self * rhs;
    }
}

impl MulAssign<&F128> for F128 {
    fn mul_assign(&mut self, rhs: &F128) {
        *self = *self * rhs;
    }
}

// Computes \sum_j COBASIS[i]^{2^j} twists[j] 
pub fn pi(i: usize, twists: &[F128]) -> F128 {
    assert!(twists.len() == 128);
    let mut ret = F128::zero();
    for j in 0..128 {
        ret += F128::from_raw(COBASIS_FROBENIUS[j][i]) * twists[j];
    }
    ret
}

#[cfg(test)]
mod tests {
    use std::{fs::File, io::Write, path::Path};
    use rand::rngs::OsRng;

    use crate::{precompute::cobasis_table::COBASIS, utils::{Matrix, _u128_from_bits}};

    use super::*;
    #[test]
    // fn test_precomputed_frobenius_matches_direct() {
    //     let diag = Matrix::diag();
    //     for k in 0..128 {
    //         for j in 0..128 {
    //             let basis_j = F128::from_binius_field(diag.cols[j]);
    //             let expected = basis_j.frob_direct(k); // Implement frob_direct as above
    //             let precomputed = FROBENIUS[k as usize][j];
    //             assert_eq!(precomputed, expected.inner_binius_field.val(), "Precomputed Frobenius mismatch at k={}, j={}", k, j);
    //         }
    //     }
    // }
    // fn test_frob_direct() {
    //     let diag = Matrix::diag();
    //     for j in 0..128 {
    //         let basis_j = F128::from_binius_field(diag.cols[j]);
    //         let mut expected = basis_j;
    //         for k in 0..128 {
    //             expected = expected.square();
    //             let direct_frob = basis_j.frob_direct(k);
    //             assert_eq!(direct_frob, expected, "frob_direct failed for k={}, j={}", k, j);
    //         }
    //     }
    // }
    #[test]
    fn test_bit_ordering() {
        let diag = Matrix::diag();
        let basis_j = F128::from_binius_field(diag.cols[0]); // Should have only bit 0 set
    
        let vec_bits = binary128_to_bits(basis_j.inner_binius_field);
        assert!(vec_bits[0], "Bit 0 should be set for the first basis vector");
        for i in 1..128 {
            assert!(!vec_bits[i], "Bit {} should not be set for the first basis vector", i);
        }
    }


    #[test]
    fn precompute_frobenius() {
        let mut basis = Matrix::diag();
        let mut ret = Vec::with_capacity(128);
        for i in 0..128 {
            ret.push(basis.cols.clone());
            for j in 0..128 {

                let x = F128::from_binius_field(basis.cols[j]);
                basis.cols[j] = x.inner_binius_field * x.inner_binius_field;
                println!("basis.cols[j] {:?}", basis.cols[j].val());
                println!("FROBENIUS[i][j] {:?}", FROBENIUS[i][j]);
                // assert_eq!(basis.cols[j].to_underlier(),FROBENIUS[i][j]);
            }   

            // // check the data just made
            // for j in 0..128 {
            //     let x = F128::from_binius_field(basis.cols[j]);
            //     let y = F128::from_binius_field(Matrix::diag().cols[j]);
            //     assert_eq!(x, y);
            // }
        }
        assert_eq!(basis, Matrix::diag());
        let ret: Vec<Vec<u128>> = ret.iter()
            .map(|row| row.iter()
            .map(|val| val.val())
            .collect())
            .collect();

        let path = Path::new("frobenius_table.txt");
        if path.is_file() {return};
        let mut file = File::create(path).unwrap();

        file.write_all("pub const FROBENIUS : [[u128; 128]; 128] =\n".as_bytes()).unwrap();
        file.write_all(format!("{:?}", ret).as_bytes()).unwrap();
        file.write_all(";".as_bytes()).unwrap();
    }

    #[test]
    fn precompute_cobasis_frobenius() {
        let path: &Path = Path::new("cobasis_frobenius_table.txt");
        if path.is_file() {return};
        let mut file = File::create(path).unwrap();
        let cobasis_vec: Vec<BinaryField128b> = COBASIS
            .iter()
            .map(|&val| BinaryField128b::new(val))
            .collect();
        let mut cobasis = Matrix::new(cobasis_vec);
        let mut ret = Vec::with_capacity(128);
        for _ in 0..128 {
            ret.push(cobasis.cols.clone());
            for j in 0..128 {
                let x = F128::from_binius_field(cobasis.cols[j]);
                cobasis.cols[j] = x.inner_binius_field * x.inner_binius_field;
            }
        }
        let ret: Vec<Vec<u128>> = ret.iter()
            .map(|row| row.iter()
            .map(|val| val.val())
            .collect())
            .collect();

        file.write_all("pub const COBASIS_FROBENIUS : [[u128; 128]; 128] =\n".as_bytes()).unwrap();
        file.write_all(format!("{:?}", ret).as_bytes()).unwrap();
        file.write_all(";".as_bytes()).unwrap();

    }

    #[test]
    fn precompute_cobasis() {
        let path = Path::new("cobasis_table.txt");
        // if path.is_file() {return};
        let mut file = File::create(path).unwrap();
        let mut matrix = vec![vec![false; 128]; 128];        
        for i in 0..128 {
            // compute pi_i linear function
            for j in 0..128 {
                let b_j = F128::basis(j);
                let b_i = F128::basis(i);
                let mut x = b_j * b_i;

                let mut s = F128::zero();
                for k in 0..128 {
                    s += x;
                    x *= x;
                }

                if s == F128::zero() {
                } else if s == F128::one() {
                        matrix[i][j] = true;
                } else {panic!()}

            }
        }
        let matrix: Vec<BinaryField128b> = matrix.iter().map(|v| BinaryField128b::new( _u128_from_bits(v))).collect();
        // let binius_matrix = binius_math::Matrix::new(matrix);
     //   let binius_matrix = binius_math::Matrix::new(matrix);
        let matrix = Matrix::new(matrix);
        let ret = matrix.inverse().unwrap().cols;
        let ret: Vec<u128> = ret.iter().map(|x| x.val()).collect::<Vec<u128>>();
        file.write_all("pub const COBASIS : [u128; 128] =\n".as_bytes()).unwrap();
        file.write_all(format!("{:?}", ret).as_bytes()).unwrap();
        file.write_all(";".as_bytes()).unwrap();
    }

    #[test]
    fn f128_is_field() {
        let rng = &mut OsRng;
        let a = F128::rand(rng);
        let b = F128::rand(rng);
        let c = F128::rand(rng);

        let one = F128::one();

        assert_eq!(a * one, a);

        assert_eq!(a + b, b + a);
        assert_eq!(a * b, b * a);
        assert_eq!((a + b) * c, a * c + b * c);
        assert_eq!((a * b) * c, a * (b * c));

        let fr = |x: F128| {x * x};

        assert_eq!(fr(a) + fr(b), fr(a + b));

        let mut x = a;
        for _ in 0..128 {
            x = fr(x);    
        }
        assert_eq!(a, x);
    }

    #[test]
    fn frobenius() {
        // pick just column i from diag
        // let col_i = F128::from_binius_field( diag.cols[i] );
        // Now check col_i.frob(k) == col_i^(2^k)

        // let rng = &mut OsRng;
        // let a = F128::rand(rng);
        let diagnal = Matrix::diag();
        let x = diagnal.cols[0];
        let y = F128::from_binius_field(x);

        assert_eq!(y, F128::one());

        for i in 0..128 {
            let col_i = F128::from_binius_field(diagnal.cols[i]);
            let mut apow = col_i;     // start apow = col_i
            for j in 0..128 {
                assert_eq!(col_i.frob(j), apow);
                apow *= apow;         // apow = apow^2, so after j steps => col_i^(2^j)
            }
        }
        // let mut apow = a;
        // for i in 0..128 {
        //     assert_eq!(a.frob(i), apow);
        //     apow *= apow;
        // }
    }

    #[test]
    fn pi_as_expected() {
        let rng = &mut OsRng;
        let mut r = F128::rand(rng);
        let mut orbit = vec![];
        for _ in 0..128 {
            orbit.push(r);
            r *= r;
        }

        for i in 0..128 {
            let lhs;
            let bit = pi(i, &orbit);
            if bit == F128::zero() {
                lhs = 0;
            } else if bit == F128::one() {
                lhs = 1;
            } else {
                panic!();
            }
            let rhs = (r.raw() >> i) % 2;
            assert!(lhs == rhs);
        }
    }

    #[test]
    fn twists_logic_and() {
        let rng = &mut OsRng;
        let a = F128::rand(rng);
        let b = F128::rand(rng);
        let mut _a = a;
        let mut _b = b;
        let mut a_orbit = vec![];
        let mut b_orbit = vec![];
        for _ in 0..128 {
            a_orbit.push(_a);
            b_orbit.push(_b);
            _a *= _a;
            _b *= _b;
        }
        let mut answer = F128::zero();
        for i in 0..128 {
            let pi_i_a = pi(i, &a_orbit);
            let pi_i_b = pi(i, &b_orbit);
            answer += F128::basis(i) * pi_i_a * pi_i_b;
        }
        let expected_answer = a & b;
        assert_eq!(answer, expected_answer);
    }
}