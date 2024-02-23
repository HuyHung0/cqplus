//! Implement structure for a table containing: size, values, and value_index_mapping
use std::collections::BTreeMap;

use ark_ff::{FftField, ToBytes};

use crate::error::Error;

#[derive(Debug)]
pub struct Table<F: FftField> {
    pub size: usize,
    pub values: Vec<F>,
    pub value_index_mapping: BTreeMap<F, usize>,
}

impl<F: FftField> ToBytes for Table<F> {
    fn write<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        self.values.write(&mut w)
    }
}

impl<F: FftField> Table<F> {
    pub fn new(values: &Vec<F>) -> Result<Self, Error> {
        if !values.len().is_power_of_two() {
            return Err(Error::TableSizeNotPow2(values.len()));
        }
        let mut value_index_mapping = BTreeMap::<F, usize>::default();
        for (i, &ti) in values.iter().enumerate() {
            let prev = value_index_mapping.insert(ti, i);
            if prev.is_some() {
                return Err(Error::DuplicateValueInTable(format!("{}", ti)));
            }
        }
        Ok(Self {
            size: values.len(),
            values: values.clone(),
            value_index_mapping,
        })
    }
}
