use crate::word::*;
use crate::errors::*;

pub fn randomize(word: &mut Word, rng: &mut dyn SeedableRng) {
    for i in 0..8 {
        let mut bytes = [0u8; 4];
        rng.fetch(&mut bytes[..]);
        word[i]=u32::from_le_bytes(bytes);
    }
}

pub fn bounded_random(top: &Word, rng: &mut dyn SeedableRng) -> Result<Word,KeygenError> {
    let mut retries = 0;
    let mut result = Word::new();
    let num_bits = top.num_bits();
    let num_words = top.num_digits();

    loop {
        randomize(&mut result, rng);
        // Mask to max have the same number of bits as top
        for i in num_words..8 {
            result[i] = 0;
        }
        result[num_words-1] &= 0xffffffff >> (num_words*32-num_bits);
        if !result.is_zero() && top.compare(&result) == 1 {
            return Ok(result);
        }else{
            retries += 1;
            if retries == 3 {
                return Err(KeygenError::FailedToGetProperRandomNumber);
            }
        }
    }
}

pub trait SeedableRng {
    fn seed(&mut self, data: &[u8]);
    fn fetch(&mut self, output: &mut [u8]);
}

pub trait Seed {
    fn fetch(&mut self) -> Option<u8>;
}
