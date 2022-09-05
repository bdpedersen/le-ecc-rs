

use crate::{word::*};


#[derive(Debug)]
pub struct Point {
    pub x: Word,
    pub y: Word,
}


impl Point {
    pub const fn new_with_contents(x: Word, y: Word) -> Point {
        Point{ x: x, y: y, }
    }

    pub fn new() -> Point {
        Point::new_with_contents(Word::new(), Word::new())
    }

    pub fn load_raw_le(&mut self, bytes: &[u8]){
        let len = bytes.len()/2;
        self.x = Word::new_from_le(&bytes[0..len]);
        self.y = Word::new_from_le(&bytes[len..2*len]);
    }

    pub fn get_raw_le(&mut self, bytes: &mut [u8]){
        let len = bytes.len()/2;
        self.x.get_raw_le(&mut bytes[0..len]);
        self.y.get_raw_le(&mut bytes[len..2*len]);
    }

    pub fn is_zero(&self) -> bool {
        self.x.is_zero() & self.y.is_zero()
    }
    // Math operations
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl Eq for Point {}

