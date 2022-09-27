// Copyright 2022, Brian Dam Pedersen.
// Copyright 2014, Kenneth MacKay
// Licensed under the BSD 2-clause license. 


use crate::word::*;
use crate::point::*;
use crate::errors::*;
use crate::rng::*;



pub trait CurveBase {

    fn mmod_fast(&self, product: &DWord) -> Word;
    fn dom_p(&self) -> &Word;
    fn dom_n(&self) -> &Word;
    fn dom_b(&self) -> &Word;
    fn dom_g(&self) -> &Point;
    fn num_bits(&self) -> u32;
}

fn bits_to_words(bits: u32) -> u32 {
    (bits+31) >> 5
}

pub trait Curve : CurveBase {
    fn new() -> Self where Self: Sized;

    fn mod_mult(&self, left: &Word, right: &Word) -> Word {
        let prod = left.mult(right);
        self.mmod_fast(&prod)
    }

    fn mod_square(&self, op: &Word) -> Word {
        let prod = op.mult(op);
        self.mmod_fast(&prod)
    }

    fn mod_add(&self, left: &Word, right: &Word) -> Word {
        left.add_mod(right, self.dom_p())
    }

    fn mod_sub(&self, left: &Word, right: &Word) -> Word {
        left.sub_mod(right, self.dom_p())
    }

    fn mod_inv(&self, op: &Word) -> Word {
        op.modular_inverse(self.dom_p())
    }

    fn double_jacobian(&self, p: &mut Point, z1: &mut Word) {
        let (x1,y1) = (&mut p.x,&mut p.y);

        let mut t4 = self.mod_square(y1); // t4 = Y1^2
        let mut t5 = self.mod_mult(x1, &t4); // t5=x1*y1^2 = A
        t4 = self.mod_square(&t4); // t4=y1^4
        *y1 = self.mod_mult(y1, z1); // t2 = y1*z1 = z3
        *z1 = self.mod_square(z1); // t3 = z1^2

        *x1 = self.mod_add(x1, z1); // t1 = x1 + z1^2
        *z1 = self.mod_add(z1,z1); // t3 = 2*z1^2
        *z1 = self.mod_sub(x1,z1); // t3 = x1-z1^2
        *x1 = self.mod_mult(x1,z1); // t1=x^2-z1^4 

        *z1 = self.mod_add(x1,x1); // t3 =2*(x1^2-z1^4)
        *x1 = self.mod_add(x1,z1); // t1 = 3*(x1^2-z14)
        // Do arithmetic shift right if x1 is odd. This is a different expression of that than the original code
        if x1.is_even() {
            x1.rshift1_inline();
        } else {
            let res = x1.add_carry(self.dom_p());
            *x1 = res.0;
            x1.rshift1_inline();
            x1[7] |= (res.1 as u32) << 31;
        }

        *z1 = self.mod_square(x1); // t3 = B^2
        *z1 = self.mod_sub(z1, &t5); // t5 = B^2 - A
        *z1 = self.mod_sub(z1, &t5); // t5 = B^2 - 2A
        t5 =self.mod_sub(&t5, z1); // t5 = A - x3
        *x1 = self.mod_mult(x1, &t5); // t1 = B* (A-x3)
        t4 = self.mod_sub(x1, &t4); // t4 = B * (A - x3) - y1^4 = y3

        *x1 = *z1;
        *z1 = *y1;
        *y1 = t4;


    }

    fn x_side(&self, x: &Word) -> Word {
        let number_three = [3u32];
        let three = Word::new_from_slice(&number_three);

        let mut result = self.mod_square(x);
        result = self.mod_sub(&result, &three);
        result = self.mod_mult(&result, x);
        self.mod_add(&result,self.dom_b())
    } 

    fn contains(&self, point: &Point) -> bool {
       if point.is_zero() {
        return false;
       } 
       
       let p: &Word = self.dom_p();

       if p.cmp_unsecure(&point.x) != 1 || p.cmp_unsecure(&point.y) != 1 {
        return false;
       }

       let tmp1 = self.mod_square(&point.y);
       let tmp2 = self.x_side(&point.x);

       tmp1.equals(&tmp2)
    }

    fn regularize_k(&self, k: &Word, kr: &mut [Word;2]) -> usize {
        let num_bits = self.num_bits();
        let num_words = bits_to_words(num_bits);
        let mut carry;
        (kr[0],carry) = k.add_carry(self.dom_n());
        carry |= (num_bits < num_words * 32) && kr[0].test_bit(num_bits as usize);
        kr[1] = kr[0].add(self.dom_n());
        
        carry as usize
    }

    // Multiplication helper functions on curve
    fn xycz_initial_double(&self, points: &mut [Point; 2],  initial_z: Option<&Word>) {
        let mut iter = points.iter_mut();
        let point2 = iter.next().unwrap();
        let point1 = iter.next().unwrap();

        let mut z: Word;
        match initial_z {
            Some(val) => {
                z = Word::new();
                z.clone_from(val);
            },
            None => {        
                let one_val = [1u32];
                z=Word::new_from_slice(&one_val[..]);
            },
        }
        point2.x = point1.x;
        point2.y = point1.y;
        self.apply_z( point1, &mut z);
        self.double_jacobian(point1, &mut z);
        self.apply_z(point2, &mut z);
    }

    fn xycz_add(&self, points: &mut [Point; 2], nb: usize) {

        // Necessary gymnastics to do branch-less swap of data
        let point_arr = [&mut points[0] as *mut Point, &mut points[1] as *mut Point]; 
        
        // Unsafe as pointers may not have been initialized - but they are provably initialized above
        let point1 = unsafe { point_arr[nb].as_mut().unwrap_unchecked() };
        let point2 = unsafe { point_arr[1-nb].as_mut().unwrap_unchecked() };
    
        let (x1,y1) = (&mut point1.x, &mut point1.y);
        let (x2,y2) = (&mut point2.x, &mut point2.y);


        let mut t5 = self.mod_sub(x2,x1);  /* t5 = x2 - x1 */
        t5 = self.mod_square(&t5);  /* t5 = (x2 - x1)^2 = A */
        *x1 = self.mod_mult(x1, &t5);              /* t1 = x1*A = B */
        *x2 = self.mod_mult(x2,&t5);      /* t3 = x2*A = C */
        *y2 = self.mod_sub(y2,y1); /* t4 = y2 - y1 */
        t5 = self.mod_square(y2);    /* t5 = (y2 - y1)^2 = D */
    
        t5 = self.mod_sub(&t5,x1); /* t5 = D - B */
        t5 = self.mod_sub(&t5,x2); /* t5 = D - B - C = x3 */
        *x2 = self.mod_sub(x2,x1); /* t3 = C - B */
        *y1 = self.mod_mult(y1,x2);   /* t2 = y1*(C - B) */
        *x2 = self.mod_sub(x1,&t5);  /* t3 = B - x3 */
        *y2 = self.mod_mult(y2, x2);               /* t4 = (y2 - y1)*(B - x3) */
        *y2 = self.mod_sub(y2,y1); /* t4 = y3 */
    
        *x2 = t5;
    }


    fn xycz_addc(&self, points: &mut [Point; 2], nb: usize) {
        // Necessary gymnastics to do branch-less swap of data
        let point_arr = [&mut points[0] as *mut Point, &mut points[1] as *mut Point]; 
        
        // Unsafe as pointers may not have been initialized - but they are provably initialized above
        let point1 = unsafe { point_arr[1-nb].as_mut().unwrap_unchecked() };
        let point2 = unsafe { point_arr[nb].as_mut().unwrap_unchecked() };
    
        let (x1,y1) = (&mut point1.x, &mut point1.y);
        let (x2,y2) = (&mut point2.x, &mut point2.y);
    
        let mut t5 = self.mod_sub(x2,x1); /* t5 = x2 - x1 */
        t5 = self.mod_square(&t5);    /* t5 = (x2 - x1)^2 = A */
        *x1 = self.mod_mult(x1, &t5);     /* t1 = x1*A = B */
        *x2 = self.mod_mult(x2, &t5);    /* t3 = x2*A = C */
        t5 = self.mod_add(y2, y1); /* t5 = y2 + y1 */
        *y2 = self.mod_sub(y2, y1); /* t4 = y2 - y1 */
    
        let mut t6 = self.mod_sub(x2, x1); /* t6 = C - B */
        *y1 = self.mod_mult(y1, &t6);   /* t2 = y1 * (C - B) = E */
        t6 = self.mod_add(x1, x2); /* t6 = B + C */
        *x2 = self.mod_square(y2);    /* t3 = (y2 - y1)^2 = D */
        *x2 = self.mod_sub(x2, &t6); /* t3 = D - (B + C) = x3 */
    
        let mut t7 = self.mod_sub(x1, x2); /* t7 = B - x3 */
        *y2 = self.mod_mult(y2, &t7);          /* t4 = (y2 - y1)*(B - x3) */
        *y2 = self.mod_sub(y2, y1); /* t4 = (y2 - y1)*(B - x3) - E = y3 */
    
        t7 = self.mod_square(&t5);   /* t7 = (y2 + y1)^2 = F */
        t7 = self.mod_sub(&t7, &t6); /* t7 = F - (B + C) = x3' */
        t6 = self.mod_sub(&t7, x1); /* t6 = x3' - B */
        t6 = self.mod_mult(&t6, &t5);      /* t6 = (y2+y1)*(x3' - B) */
        *y1 = self.mod_sub(&t6, y1); /* t2 = (y2+y1)*(x3' - B) - E = y3' */
    
        *x1 = t7; 
    

    }

    fn apply_z(&self, point: &mut Point, z: &Word) {
        let (x1,y1) = (&mut point.x, &mut point.y);

        let mut t1 = self.mod_square(z);
        *x1 = self.mod_mult(x1, &t1);
        t1 = self.mod_mult(&t1, &z);
        *y1 = self.mod_mult(y1, &t1);
    }

    fn mult(&self, point: &Point, scalar: &Word, initial_z: Option<&Word>) -> Point {
        let mut r = [Point::new(), Point::new()];
        r[1].x = point.x;
        r[1].y = point.y;

        self.xycz_initial_double(&mut r, initial_z);
        // NOTE: original code implicitly adds 1 bit to num_bits here !!!
        for i in (1..self.num_bits()).rev() {
            let nb = (!scalar.test_bit(i as usize)) as usize;
            self.xycz_addc(&mut r, nb);
            self.xycz_add(&mut r, nb);
        }

        let nb = (!scalar.test_bit(0)) as usize;
        self.xycz_addc(&mut r, nb);

        let mut z = self.mod_sub(&r[1].x, &r[0].x);
        z = self.mod_mult(&z, &r[1-nb].y);
        z = self.mod_mult(&z, &point.x);
        z = self.mod_inv(&z);
        z = self.mod_mult(&z, &point.y);
        z = self.mod_mult(&z, &r[1-nb].x);

        self.xycz_add(&mut r, nb);
        self.apply_z(&mut r[0],&z);

        let return_value: Point = Point {x: r[0].x, y: r[0].y, };
        return_value
    }

    fn mult_regularized(&self, public: &Point, private: &Word, rng: Option<&mut dyn SeedableRng>) -> Result<Point,KeygenError> {
        let mut tmp = [Word::new(), Word::new()];
        let mut initial_z: Option<&Word> = None;
        let base: &[u32; 8];

        tmp[0] = *private;

        let carry = self.regularize_k(private, &mut tmp);

        if let Some(myrng) = rng {
            let result = bounded_random(self.dom_p(), myrng)?;
            tmp[carry].copy_from_slice(&result[..]);
            initial_z = Some(&tmp[carry]);
        }

        base = &tmp[1-carry];
        // Compute the actual result
        let result = self.mult(public,base, initial_z);
        if result.is_zero() {
            Err(KeygenError::PointAtInfinity)
        } else {
            Ok(result)
        }

    }

    fn compute_public_key(&self, private: &Word, rng: Option<&mut dyn SeedableRng>) -> Result<Point,KeygenError> {
        self.mult_regularized(self.dom_g(), private, rng)
     }

     fn compute_shared_secret(&self, own_private: &Word, other_public: &Point, rng: Option<&mut dyn SeedableRng>) -> Result<Point,KeygenError> {
        self.mult_regularized(other_public, own_private, rng)
     }

}

pub struct CurveSecp256r1 {
    p: Word,
    n: Word,
    g: Point,
    b: Word, 
    num_bits: u32,
}

impl Curve for CurveSecp256r1 {
    fn new() -> Self {
        let p: [u8; 32] = [
            0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
            0xff,0xff,0xff,0xff,0x00,0x00,0x00,0x00,
            0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
            0x01,0x00,0x00,0x00,0xff,0xff,0xff,0xff,
            ];
        let n: [u8; 32] = [
            0x51,0x25,0x63,0xfc,0xc2,0xca,0xb9,0xf3,
            0x84,0x9e,0x17,0xa7,0xad,0xfa,0xe6,0xbc,
            0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
            0x00,0x00,0x00,0x00,0xff,0xff,0xff,0xff,
            ];
        let gx: [u8; 32] = [
            // x
            0x96, 0xc2, 0x98, 0xd8, 0x45, 0x39, 0xa1, 0xf4,
            0xa0, 0x33, 0xeb, 0x2d, 0x81, 0x7d, 0x03, 0x77,
            0xf2, 0x40, 0xa4, 0x63, 0xe5, 0xe6, 0xbc, 0xf8,
            0x47, 0x42, 0x2c, 0xe1, 0xf2, 0xd1, 0x17, 0x6b,
            ];
            // y
        let gy: [u8; 32] = [
            0xf5, 0x51, 0xbf, 0x37, 0x68, 0x40, 0xb6, 0xcb,
            0xce, 0x5e, 0x31, 0x6b, 0x57, 0x33, 0xce, 0x2b,
            0x16, 0x9e, 0x0f, 0x7c, 0x4a, 0xeb, 0xe7, 0x8e,
            0x9b, 0x7f, 0x1a, 0xfe, 0xe2, 0x42, 0xe3, 0x4f,
            ];
        let b: [u8; 32] = [
            0x4b, 0x60, 0xd2, 0x27, 0x3e, 0x3c, 0xce, 0x3b,
            0xf6, 0xb0, 0x53, 0xcc, 0xb0, 0x06, 0x1d, 0x65,
            0xbc, 0x86, 0x98, 0x76, 0x55, 0xbd, 0xeb, 0xb3,
            0xe7, 0x93, 0x3a, 0xaa, 0xd8, 0x35, 0xc6, 0x5a,            
            ];

        CurveSecp256r1 { 
            p: Word::new_from_le(&p), 
            n: Word::new_from_le(&n), 
            g: Point::new_with_contents(Word::new_from_le(&gx), Word::new_from_le(&gy)), 
            b: Word::new_from_le(&b), 
            num_bits: 256 
        }
    }
}

impl CurveBase for CurveSecp256r1 {


    fn dom_p(&self) -> &Word {
        &self.p
    }

    fn dom_b(&self) -> &Word {
        &self.b
    }

    fn dom_g(&self) -> &Point {
        &self.g
    }

    fn dom_n(&self) -> &Word {
        &self.n
    }

    fn num_bits(&self) -> u32 {
        self.num_bits
    }
    
    fn mmod_fast(&self, product: &DWord) -> Word {
        // Note: in the original code we substract carry from sub-methods, but my sub-methods return
        // -1 for borrow so they need to add

        let mut result = Word::new_from_slice(&product[0..8]);
        let mut tmp = Word::new();

        // S1
        for i in 3..8 {
            tmp[i]=product[i+8];
        }

        let mut carry = tmp.double_inline();
        carry += result.add_inline(&tmp);

        // S2
        for i in 3..7 {
            tmp[i]=product[i+9];
        }
        tmp[7]=0;

        carry += tmp.double_inline();
        carry += result.add_inline(&tmp);

        // S3
        for i in 0..8 {
            tmp[i]=product[i+8];
        }

        for i in 3..6 {
            tmp[i] = 0;
        }

        carry += result.add_inline(&tmp);

        // s4

        for i in 0..3 {
            tmp[i] = product[i+9];
        }

        for i in 3..6 {
            tmp[i] = product[i+10];
        }

        tmp[6] = product[13];
        tmp[7] = product[8];
        carry += result.add_inline(&tmp);

        // d1
        for i in 0..3 {
            tmp[i]=product[i+11];
        }

        for i in 3..6 {
            tmp[i] = 0;
        }

        tmp[6] = product[8];
        tmp[7] = product[10];
        carry += result.sub_inline(&tmp);

        // d2
        for i in 0..4 {
            tmp[i]=product[i+12];
        }

        for i in 4..6 {
            tmp[i] = 0;
        }

        tmp[6] = product[9];
        tmp[7] = product[11];
        carry += result.sub_inline(&tmp);

        // d3
        for i in 0..3 {
            tmp[i]=product[i+13];
        }

        for i in 3..6 {
            tmp[i] = product[i+5];
        }

        tmp[6] = 0;
        tmp[7] = product[12];
        carry += result.sub_inline(&tmp);

        // d4
        for i in 0..2 {
            tmp[i]=product[i+14];
        }
        tmp[2] = 0;

        for i in 3..6 {
            tmp[i] = product[i+6];
        }

        tmp[6] = 0;
        tmp[7] = product[13];
        carry += result.sub_inline(&tmp);

        if carry < 0 {
            loop {
                carry += result.add_inline(&self.p);
                if carry == 0 {break;}
            }
        }else {
            while carry != 0 || self.p.cmp_unsecure(&result) != 1 {
                carry += result.sub_inline(&self.p);
            }
        }

        result
    } 
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::le_testdata::*;

    struct DilbertRng {}

    impl SeedableRng for DilbertRng {
        fn seed(&mut self, _data: &[u8]) {}
        fn fetch(&mut self, output: &mut [u8]) {
            for i in 0..output.len() {
                output[i]=9;
            }
        }
    }

    #[test]
    fn check_le_data_pub_a() -> Result<(),KeygenError> {
        let curve = CurveSecp256r1::new();
        // Check public key A
        for i in 0..2 {
            assert!(curve.contains(&le_testdata::PUB_A[i]));
        }
        Ok(())
    }

    #[test]
    fn check_le_data_pub_b() -> Result<(),KeygenError> {
        let curve = CurveSecp256r1::new();
        // Check public key B
        for i in 0..2 {
            assert!(curve.contains(&le_testdata::PUB_B[i]));
        }
        Ok(())
    }

    #[test]
    fn check_le_data_priv_a() -> Result<(),KeygenError> {
        let curve = CurveSecp256r1::new();
        // Check that the private keys in A generates the public keys
        for i in 0..2 {
            let pubkey = curve.compute_public_key(&le_testdata::PRIV_A[i], None)?;
            assert_eq!(pubkey,le_testdata::PUB_A[i]);
        }
        Ok(())
    }

    #[test]
    fn check_le_data_priv_b() -> Result<(),KeygenError> {
        let curve = CurveSecp256r1::new();
        // Check that the private keys in B generates the public keys
        for i in 0..2 {
            let pubkey = curve.compute_public_key(&le_testdata::PRIV_B[i], None)?;
            assert_eq!(pubkey,le_testdata::PUB_B[i]);
        }

        Ok(())
    }

    #[test]
    fn check_le_data_dh_ab() -> Result<(), KeygenError> {
        let curve = CurveSecp256r1::new();
        for i in 0..2 {
            let shs = curve.compute_shared_secret(&le_testdata::PRIV_A[i], 
                &le_testdata::PUB_B[i], None)?;
            assert_eq!(shs.x,le_testdata::DHKEY[i]);
        }

        Ok(())
    }

    #[test]
    fn check_le_data_dh_ba() -> Result<(), KeygenError> {
        let curve = CurveSecp256r1::new();
        for i in 0..2 {
            let shs = curve.compute_shared_secret(&le_testdata::PRIV_B[i], 
                &le_testdata::PUB_A[i], None)?;
            assert_eq!(shs.x,le_testdata::DHKEY[i]);
        }

        Ok(())
    }

    #[test]
    fn check_le_data_dh_ab_rng() -> Result<(), KeygenError> {
        let curve = CurveSecp256r1::new();
        let mut rng = DilbertRng {};
        for i in 0..2 {
            let shs = curve.compute_shared_secret(&le_testdata::PRIV_A[i], 
                &le_testdata::PUB_B[i], Some(&mut rng))?;
            assert_eq!(shs.x,le_testdata::DHKEY[i]);
        }

        Ok(())
    }

    #[test]
    fn check_le_data_dh_ba_rng() -> Result<(), KeygenError> {
        let curve = CurveSecp256r1::new();
        let mut rng = DilbertRng {};
        for i in 0..2 {
            let shs = curve.compute_shared_secret(&le_testdata::PRIV_B[i], 
                &le_testdata::PUB_A[i], Some(&mut rng))?;
            assert_eq!(shs.x,le_testdata::DHKEY[i]);
        }

        Ok(())
    }


}

