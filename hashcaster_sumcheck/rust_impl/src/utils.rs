// This mainly contains linear algebra code required to find necessary bases.

use binius_field::{underlier::WithUnderlier, BinaryField128b};
use bytemuck::cast;
use rand::Rng;

use crate::field::F128;

pub fn log2_exact(mut x: usize) -> usize {
    let mut c = 0;
    loop {
        if x % 2 == 1 {
            if x == 1 {
                return c;
            } else {
                panic!("Not an exact power of 2.");
            }
        }
        c += 1;
        x >>= 1;
    }
}

pub fn binary128_to_bits(x: BinaryField128b) -> Vec<bool> {
    let mut ret = Vec::with_capacity(128);
    let bytes = cast::<u128, [u8; 16]>(x.val());
    
    for i in 0..16 {
        for j in 0..8 {
            ret.push((bytes[i] & (1 << j)) != 0)
        }
    }
    ret
}

pub fn u128_idx(x: &u128, i: usize) -> bool {
    let bytes = cast::<u128, [u8; 16]>(*x);
    let byte_idx = i >> 3;
    let bit_idx = i ^ (byte_idx << 3);
    bytes[byte_idx] & (1 << bit_idx) != 0
}

/// Not efficient.
pub fn _u128_from_bits(bits: &[bool]) -> u128 {
    let l = bits.len();
    let mut ret = 0;
    assert!(l <= 128);
    for i in 0..l {
        ret += (bits[i] as u128) << i
    }
    ret
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix {
    pub cols: Vec<BinaryField128b>,
}

impl Matrix {
    pub fn new(cols: Vec<BinaryField128b>) -> Self {
        assert_eq!(cols.len(), 128);
        Self{ cols }
    }

/// Applies 128x128 matrix in column form to 128 vector.
/// Somewhat efficient.
    pub fn apply(&self, vec: BinaryField128b) -> BinaryField128b {
        let matrix: &Vec<BinaryField128b> = &self.cols;
        let mut ret = BinaryField128b::new(0);
        let vec_bits = binary128_to_bits(vec);
        for i in 0..128 {
            if vec_bits[i] {ret += matrix[i]}
        }
        ret
    }

    pub fn compose(&self, b: &Self) -> Self {
        let b = &b.cols;
        let mut ret = Vec::with_capacity(128);
        for i in 0..128 {
            ret.push(self.apply(b[i]));
        }
        Self{cols: ret}
    }

    pub fn diag() -> Self {
        let mut ret = Vec::with_capacity(128);
        let mut a = 1;
        for i in 0..128 {
            let item = F128::basis(i).inner_binius_field;
            // assert!(a == item.val());
            ret.push(item);
            a <<= 1; 

            // println!("{:b}", item.val());
        }
        Self::new(ret)
    }

    pub fn swap_cols(&mut self, i: usize, j: usize) {
        let tmp = self.cols[j];
        self.cols[j] = self.cols[i];
        self.cols[i] = tmp;
    }

    pub fn triang(&mut self, i: usize, j: usize) {
        let temp = self.cols[i];
        self.cols[j] += temp;
    }

    pub fn inverse(&self) -> Option<Self> {
        let mut a = self.clone();
        let mut b = Self::diag();
        for i in 0..128 {
            let mut j = i;
            if !u128_idx(&a.cols[i].val(), i) {
                loop {
                    j += 1;
                    if j == 128 {
                        return None;
                    }
                    if u128_idx(&a.cols[j].val(), i) {
                        a.swap_cols(i, j);
                        b.swap_cols(i, j);
                        break;
                    }

                }
            }
            loop {
                j += 1;
                if j == 128 {break}
                if u128_idx(&a.cols[j].val(), i) {
                    a.triang(i, j);
                    b.triang(i, j);
                }
            }
        }
        // a is now under-diagonal
        for i in 1..128 {
            for j in 0..i {
                if u128_idx(&a.cols[j].val(), i) {
                    a.triang(i, j);
                    b.triang(i, j);
                }
            }
        }
        Some(b)
    }
    
}

pub fn u128_rand<RNG: Rng>(rng: &mut RNG) -> u128 {
    let a = rng.next_u64();
    let b = rng.next_u64();
    cast::<[u64; 2], u128>([a, b])
}

#[cfg(test)]
mod tests {
    use rand::{rngs::OsRng, RngCore};
    use super::*;

    #[test]
    fn test_apply_matrix() {
        let rng = &mut OsRng;
        let mut matrix = vec![];
        for _ in 0..128 {
            matrix.push( BinaryField128b::new( u128_rand(rng)));
        }
        let vec = u128_rand(rng);
        let vec_bits = binary128_to_bits(BinaryField128b::new( vec));
        let matrix_bits : Vec<Vec<bool>> = matrix.iter().map(|v|binary128_to_bits(*v)).collect();
        let mut expected_answer = vec![false; 128];
        for i in 0..128 {
            for j in 0..128 {
                expected_answer[j] ^= matrix_bits[i][j] && vec_bits[i]
            }
        }
        let answer = Matrix::new(matrix).apply(BinaryField128b::new((vec)));
        assert_eq!(binary128_to_bits(answer), expected_answer);
    }
    #[test]
    fn test_invert_matrix() {
        let rng = &mut OsRng;
        let mut matrix = Matrix::diag();
        // generate (kind of) random invertible matrix by applying a lot of triangular transforms and swaps
        for _ in 0..100_000 {
            let r = rng.next_u64() as usize;
            let r1 = r % (1 << 4);
            let r2 = (r >> 4) & 127;
            let r3 = (r >> 32) & 127;
            if r1 == 0 {
                matrix.swap_cols(r2, r3);
            } else {
                if r2 != r3 {
                    matrix.triang(r2, r3);
                }
            }
        }
        let test_vector = BinaryField128b::new(u128_rand(rng));
        let inv = matrix.inverse().unwrap();
        assert_eq!(
            matrix.apply(inv.apply( test_vector )),
            test_vector,
        );
        assert_eq!(
            inv.compose(&matrix),
            Matrix::diag(),
        )
    }

}