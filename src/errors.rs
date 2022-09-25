// Copyright 2022, Brian Dam Pedersen.
// Copyright 2014, Kenneth MacKay
// Licensed under the BSD 2-clause license. 

use core::fmt;

#[derive(Debug)]
pub enum KeygenError {
    NotEnoughEntropy,
    PointAtInfinity,
    FailedToGetProperRandomNumber,
    KeyNotOnCurve,
    NoPrivateKeyAvailable,
}

impl fmt::Display for KeygenError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            KeygenError::NotEnoughEntropy => write!(f, "Not enough Entropy"),
            KeygenError::PointAtInfinity => write!(f, "Point at infinity"),
            KeygenError::FailedToGetProperRandomNumber => write!(f, "Failed to get proper random number"),
            KeygenError::KeyNotOnCurve => write!(f,"Key not on curve"),
            KeygenError::NoPrivateKeyAvailable => write!(f,"No private key available"),
        }
    }

}
