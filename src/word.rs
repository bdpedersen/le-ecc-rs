//pub mod word {
use core::slice::{IterMut};
use core::assert;


pub type Word = [u32; 8];
pub type DWord = [u32; 16];

pub trait WordConstruction {
    fn new() -> Word; // *
    fn new_from_slice(val: &[u32]) -> Word; // *
    fn new_from_le(bytes: &[u8]) -> Word;
    fn get_raw_le(&self, bytes: &mut [u8]);
}




impl WordConstruction for Word {
    fn new() -> Word {
        [0; 8]
    }
    
    fn new_from_slice(val: &[u32]) -> Word {
        let mut r = Word::new();
        let mut iter = val.iter();
        for i in 0..8 {
            if let Some(v) = iter.next() {
                r[i] = *v;
            }else{
                break;
            }
        }
        r
    }
    
    fn new_from_le(bytes: &[u8]) -> Word {
        let mut buf = Word::new();
            assert!(bytes.len()*(u8::BITS as usize) <= 8*(u32::BITS as usize) && bytes.len() % 4 == 0);
        for i in 0..buf.len() {
            let bidx = i*4;
            let arr: [u8; 4] = bytes[bidx..bidx+4].try_into().unwrap();
            buf[i]=u32::from_le_bytes(arr);
        }
        buf
    }
    
    
    fn get_raw_le(&self, bytes: &mut [u8]){
        assert!(bytes.len()*(u8::BITS as usize) >= 8*(u32::BITS as usize) && bytes.len() % 4 == 0);
        for i in 0..self.len() {
            let idx = self.len()-i-1;
            let bidx = i*4;
            let arr:[u8 ; 4] = self[idx].to_le_bytes();
            for j in 0..4 {
                bytes[bidx+j]= arr[j];
            }
        }
        
    }
}

pub trait WordOps : WordConstruction {
    fn cmp_unsecure(&self, other: &Word) -> i32; // *
    fn rshift1_inline(&mut self); // *
    fn ari_rshift1_inline(&mut self); // *
    fn equals(&self, other: &Word) -> bool; // constant-time so secure
    fn add_inline(&mut self, other: &Word) -> i32; // *
    fn sub_inline(&mut self, other: &Word) -> i32; // *
    fn double_inline(&mut self) -> i32; // *
    fn mult(&self, other: &Word) -> DWord; // *
    fn add_mod(&self, other: &Word, modulus: &Word) -> Word; // *
    fn sub_mod(&self, other: &Word, modulus: &Word) -> Word; // *
    fn is_even(&self) -> bool;
    fn is_zero(&self) -> bool;
    fn add_carry(&self, other: &Word) -> (Word,bool);
    fn test_bit(&self, bit: usize) -> bool;
    fn add(&self, other: &Word) -> Word;
    fn compare(&self, other: &Word) -> i32;
    fn sub_borrow(&self, other: &Word) -> (Word, bool);
    fn num_bits(&self) -> usize;
    fn num_digits(&self) -> usize;
    fn modular_inverse(&self, modulus: &Word) -> Word;
    fn modinv_update_inline(&mut self, modulus: &Word);
}

// helper functions
fn rshift_iter_inline(iter_mut: IterMut<'_,u32>) {
    let mut carry = 0u32;
    for r in iter_mut.rev() {
        let x = *r;
        *r = carry | x >> 1;
        carry = x << 31;
    }
}

fn muladd(a: u32, b: u32, r2: &mut u32, r1: &mut u32, r0: &mut u32) {
    let prod: u64 = (a as u64)*(b as u64);
    let mut accum = (prod as u32 as u64) + (*r0 as u64);
    *r0 = accum as u32;
    accum >>= 32;
    accum += (*r1 as u64) + (prod >> 32);
    *r1 = accum as u32;
        *r2 += (accum >> 32) as u32;
}


impl WordOps for Word {
    fn is_zero(&self) -> bool {
        self.iter().fold(true, |acc, elem| acc & (*elem == 0))
    }

    fn cmp_unsecure(&self, other: &Word) -> i32 {
        for i in (0..8).rev() {
            let idx = i as usize;
            if self[idx] > other[idx] {
                return 1;
            }
            if self[idx] < other[idx] {
                return -1;
            }
        }
        0
    }



    fn equals(&self, other: &Word) -> bool {
        // Only one test in the end to avoid spilling information
        self.iter().zip(other.iter()).fold(0u32, |a, x| { a | (x.0.wrapping_sub(*x.1))}) == 0
    }

    fn add_inline(&mut self, other: &Word) -> i32 {
        let mut carry: u64 = 0;
        for (r,y) in self.iter_mut().zip(other.iter()) {
            let accum = (*r as u64).wrapping_add(*y as u64).wrapping_add(carry);
            carry = accum >> 32;
            *r =  accum as u32;
        }
        carry as i32
    }

    fn sub_inline(&mut self, other: &Word) -> i32 {
        let mut borrow: u64 = 0;
        for (r,y) in self.iter_mut().zip(other.iter()) {
            let accum = (*r as u64).wrapping_sub(*y as u64).wrapping_add(borrow);
            borrow = (accum >> 32) as i32 as u64;
            *r = accum as u32;
        }
        borrow  as i32
    }

    fn double_inline(&mut self) -> i32 {
        let mut carry: u64 = 0;
        for r in self.iter_mut() {
            let y = *r;
            let accum = (y as u64).wrapping_add(y as u64).wrapping_add(carry);
            carry = accum >> 32;
            *r =  accum as u32;
        }
        carry as i32
    }

    fn rshift1_inline(&mut self) {
        rshift_iter_inline(self.iter_mut());
    }

    fn ari_rshift1_inline(&mut self) {
        let sbit = self[7] & 0x80000000;
        rshift_iter_inline(self.iter_mut());
        self[7] |= sbit;
    }

    fn mult(&self, other: &Word) -> DWord {
        let mut r0 = 0u32;
        let mut r1 = 0u32;
        let mut r2 = 0u32;
        let mut result: DWord = [0; 16];
        for k in 0..8 {
            for i in 0..=k {
                muladd(self[i],other[k-i],&mut r2, &mut r1, &mut r0);
            }
            result[k]=r0;
            r0=r1;
            r1=r2;
            r2=0;
        }
        for k in 8..15 {
            for i in (k+1-8)..8 {
                muladd(self[i],other[k-i],&mut r2, &mut r1, &mut r0);
            }
            result[k]=r0;
            r0=r1;
            r1=r2;
            r2=0;
        }
        result[15]=r0;
        result
    }

    fn add_mod(&self, other: &Word, modulus: &Word) -> Word {
        let res = self.add_carry(other);
        let mut result = res.0;
        let carry = res.1;
        if carry || modulus.cmp_unsecure(&result) != 1 {
            result.sub_inline(modulus);
        }
        result
    }

    fn sub_mod(&self, other: &Word, modulus: &Word) -> Word {
        let res = self.sub_borrow(other);
        let mut result = res.0;
        let borrow = res.1;
        if borrow  {
            result.add_inline(modulus);
        }
        result
    }

    fn is_even(&self) -> bool {
        self[0] & 1 == 0
    }

    fn add_carry(&self, other: &Word) -> (Word, bool) {
        let mut carry: u64 = 0;
        let mut result = Word::new();
        for ((x,y),r) in self.iter().zip(other.iter()).zip(result.iter_mut()) {
            let accum: u64 = (*x as u64) + (*y as u64) + carry;
            carry = accum >> 32;
            *r = accum as u32;
        }
        (result, carry != 0)
    }

    fn test_bit(&self, bit: usize) -> bool {
        self[bit >> 5] & (1 << (bit & 0x1f)) > 0
    }

    fn add(&self, other: &Word) -> Word {
        self.add_carry(other).0
    }

    fn compare(&self, other: &Word) -> i32 {
        let tmp = self.sub_borrow(other);
        let neg = tmp.1 as i32;
        let nonzero = (!tmp.0.is_zero()) as i32;
        
        nonzero - 2*neg
    }
    fn sub_borrow(&self, other: &Word) -> (Word, bool) {
        let mut borrow: u64 = 0;
        let mut result = Word::new();
        for ((x,y),r) in self.iter().zip(other.iter()).zip(result.iter_mut()) {
            let accum: u64 = (*x as u64).wrapping_sub(*y as u64).wrapping_add(borrow);
            borrow = (accum >> 32) as i32 as u64; // Sign extend carry
            *r = accum as u32;
        }
        (result, borrow != 0)
    }
    fn num_digits(&self) -> usize {
        let mut i:i32= (self.len() as i32)-1;
        while i >= 0 && self[(i as usize)]==0 {
            i-=1;
        }
        (i as usize) + 1
    }      
    
    fn num_bits(&self) -> usize {
        let num_digits = self.num_digits();
        if num_digits == 0 {
            0
        }else{
            let mut digit = self[num_digits-1];
            let mut i=0;
            while digit != 0 {
                digit >>= 1;
                i+=1;
            }
            ((num_digits-1) << 5) + i
        }
    }

    fn modinv_update_inline(&mut self, modulus: &Word) {
        let mut carry = 0;
        if !self.is_even() {
            carry = ((self.add_inline(modulus) != 0) as u32) << 31;
        }
        self.rshift1_inline();
        self[7] |= carry;
    }

    fn modular_inverse(&self, modulus: &Word) -> Word {
        if self.is_zero() {
            return Word::new();
        }

        let mut a = self.clone();
        let mut b = modulus.clone();
        let mut u = Word::new();
        u[0]=1;
        let mut v = Word::new();
        let mut cmpresult = a.cmp_unsecure(&b);
        while cmpresult != 0 {
            if a.is_even() {
                a.rshift1_inline();
                u.modinv_update_inline(modulus);
            } else if b.is_even() {
                b.rshift1_inline();
                v.modinv_update_inline(modulus);
            } else if cmpresult > 0 {
                a.sub_inline(&b);
                a.rshift1_inline();
                if u.cmp_unsecure(&v) < 0 {
                    u.add_inline(modulus);
                }
                u.sub_inline(&v);
                u.modinv_update_inline(modulus);
            } else {
                b.sub_inline(&a);
                b.rshift1_inline();
                if v.cmp_unsecure(&u) < 0 {
                    v.add_inline(modulus);
                }
                v.sub_inline(&u);
                v.modinv_update_inline(modulus);
            }
            cmpresult = a.cmp_unsecure(&b);
        }
        u
    }
    

}


trait DWordMath {
    fn new() -> DWord;
    fn set(&mut self, other: &Word, offset: usize);
}

impl DWordMath for DWord {
    fn new() -> DWord {
        [0; 16]
    }

    fn set(&mut self, other: &Word, offset: usize) {
        for i in 0..8 {
            self[i+offset]=other[i];
        }
    }


}
//}
