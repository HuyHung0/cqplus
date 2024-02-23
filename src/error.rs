//! Contain some pre-defined errors
#[derive(Debug, PartialEq)]
pub enum Error {
    TableSizeNotPow2(usize),
    WitnessSizeNotPow2(usize),
    DuplicateValueInTable(String),
    ValueNotInTable(String),
    VerificationEquation1Error,
    VerificationEquation2Error,
    CommitmentError,
}
