// Copyright 2022, Brian Dam Pedersen.
// Copyright 2014, Kenneth MacKay
// Licensed under the BSD 2-clause license. 

use crate::word::*;
use crate::point::*;
use crate::ecc_curve::*;
use crate::errors::*;
use crate::rng::*;

fn max(left: i32, right: i32) -> i32 {
    if left > right {left} else {right}
}



pub struct KeySetP256 {
    private: Option<Word>,
    public: Option<Point>,
    seed_quality: i32,
    is_debug_key: bool,
    curve: CurveSecp256r1,
}

impl KeySetP256 {
    pub fn new() -> KeySetP256 {
        KeySetP256 { private: None, public: None, seed_quality: 0, is_debug_key: false, curve: CurveSecp256r1::new() }
    }
}

mod private {
    use crate::word::Word;
    use crate::point::Point;
    use crate::ecc_curve::Curve;
    pub trait KeySetAccessors {
        fn _private(&self) -> &Option<Word>;
        fn _set_private(&mut self,val: Option<Word>);
        fn _public(&self) -> &Option<Point>;
        fn _set_public(&mut self, val: Option<Point>);
        fn _seed_quality(&self) -> i32;
        fn _set_seed_quality(&mut self, val: i32);
        fn _set_debug_flag(&mut self);
        fn _clear_debug_flag(&mut self);
        fn _debug_flag(&self) -> bool;
        fn _curve(&self) -> &dyn Curve;
    }
}

pub trait KeySet : private::KeySetAccessors {

    fn debug_set_private_key(&mut self, key: &Word) {
        self._set_private(Some(key.clone()));
        if let Ok(val) = self._curve().compute_public_key(key,  None) {
            self._set_public(Some(val));
        }else{
            panic!("Invalid debug key");
        }
        self._set_debug_flag();
    }

    fn gen_keyset(&mut self, rng: &mut dyn SeedableRng, seed: &mut dyn Seed) -> Result<(),KeygenError> {
        // Seed the random number generator with at least 8 bytes
        let mut bytecnt = 0;
        let mut seeddata = [0u8; 16];

        self._clear_debug_flag();
        for i in 0..16 {
            if let Some(x) = seed.fetch() {
                seeddata[i] = x;
                bytecnt += 1;
            }else{
                break;
            }
        }
        rng.seed(&seeddata[..bytecnt]);
        let mut seed_quality = self._seed_quality();
        seed_quality += bytecnt as i32;
        self._set_seed_quality(seed_quality);
        if seed_quality < 8 {
            return Err(KeygenError::NotEnoughEntropy);
        }

        // Set bound for private key
        let mut privkey = Word::new();
        let r = self._curve().dom_n().clone();
        while privkey.is_zero() {
            privkey =  bounded_random(&r, rng)?;
        }

        let pubkey = self._curve().compute_public_key(&privkey, Some(rng))?;
        self._set_private(Some(privkey));
        self._set_public(Some(pubkey));
        return Ok(());    
    }

    fn erase_keyset(&mut self) {
        self._set_seed_quality(max(0,self._seed_quality()-8));
        self._set_private(None);
        self._set_public(None);
        self._clear_debug_flag();
    }

    fn public_key(&self) -> &Option<Point> {
        self._public()
    }

    fn shared_secret(&self, other_public_key: &Point, rng: Option<&mut dyn SeedableRng>) -> Result<Word,KeygenError> {
        if !self._curve().contains(other_public_key) {
            return Err(KeygenError::KeyNotOnCurve);
        }

        if let None = rng {
            if !self._debug_flag() {
                return Err(KeygenError::NotEnoughEntropy);
            }
        }

        if let Some(private) = self._private() {
            let result = self._curve().compute_shared_secret(private, other_public_key, rng)?;
            return Ok(result.x);
        } else {
            return Err(KeygenError::NoPrivateKeyAvailable);
        }
    }
}

impl KeySet for KeySetP256 {}

impl private::KeySetAccessors for KeySetP256 {
    fn _private(&self) -> &Option<Word> {
        &self.private
    }

    fn _set_private(&mut self, val: Option<Word>) {
        self.private = val
    }

    fn _public(&self) -> &Option<Point> {
        &self.public
    }
    fn _set_public(&mut self, val: Option<Point>) {
        self.public = val;
    }

    fn _seed_quality(&self) -> i32 {
        self.seed_quality
    }

    fn _set_seed_quality(&mut self, val: i32) {
        self.seed_quality = val;
    }

    fn _curve(& self) -> &dyn Curve {
        &self.curve
    }

    fn _set_debug_flag(&mut self) {
        self.is_debug_key = true;
    }

    fn _clear_debug_flag(&mut self) {
        self.is_debug_key = false;
    }

    fn _debug_flag(&self) -> bool {
        self.is_debug_key
    }

}


#[cfg(test)]
extern crate krang_rs;

// Tests
#[cfg(test)]
mod tests {
    use super::*;


    struct MyRng {
        rng: krang_rs::Krang,
    }

    impl MyRng {
        fn new() -> MyRng {
            MyRng {rng: krang_rs::Krang::new() }
        }
    }

    impl SeedableRng for MyRng {
        fn seed(&mut self, data: &[u8]) {
            self.rng.seed(data);
        }

        fn fetch(&mut self, output: &mut [u8]) {
            self.rng.fetch(output);
        }
    }

    struct MySeed {}

    impl Seed for MySeed {
        fn fetch(&mut self) -> Option<u8> {
            Some(9)
        }
    }

    #[test]
    fn exchange_dh() -> Result<(), KeygenError> {
        let mut rng = MyRng::new();
        let mut seed = MySeed {};
        let mut alice = KeySetP256::new();
        let mut bob = KeySetP256::new();
        
        alice.gen_keyset(&mut rng, &mut seed)?;
        bob.gen_keyset(&mut rng, &mut seed)?;

        let bobs_public_key = bob.public_key().as_ref().unwrap();
        let alices_public_key = alice.public_key().as_ref().unwrap();

        assert!(bobs_public_key != alices_public_key);

        let bobs_secret = bob.shared_secret(alices_public_key,Some(&mut rng))?;
        let alices_secret = alice.shared_secret(bobs_public_key, Some(&mut rng))?;

        assert_eq!(bobs_secret,alices_secret);

        alice.erase_keyset();

        assert!(alice.public_key().is_none());
        alice.gen_keyset(&mut rng, &mut seed)?;
        let alices_public_key = alice.public_key().as_ref().unwrap();

        let alices_secret = alice.shared_secret(bobs_public_key, Some(&mut rng))?;

        assert!(alices_secret != bobs_secret);

        let bobs_secret = bob.shared_secret(alices_public_key,Some(&mut rng))?;

        assert_eq!(bobs_secret,alices_secret);

        let mut alices_wrong_key = alices_public_key.clone();

        alices_wrong_key.x[0]=0xdeadbeef;

        let bobs_secret = bob.shared_secret(&alices_wrong_key,Some(&mut rng));

        assert!(bobs_secret.is_err());

        alices_wrong_key = alices_public_key.clone();
        alices_wrong_key.y = Word::new();

        let bobs_secret = bob.shared_secret(&alices_wrong_key,Some(&mut rng));

        assert!(bobs_secret.is_err());

        alices_wrong_key = alices_public_key.clone();
        alices_wrong_key.x = Word::new();

        let bobs_secret = bob.shared_secret(&alices_wrong_key,Some(&mut rng));

        assert!(bobs_secret.is_err());

        alices_wrong_key = Point::new();

        let bobs_secret = bob.shared_secret(&alices_wrong_key,Some(&mut rng));

        assert!(bobs_secret.is_err());



        Ok(())
    }
}
