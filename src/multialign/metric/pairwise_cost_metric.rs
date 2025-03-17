use std::{marker::PhantomData, path::Path};

use anyhow::{Context, Result};
use compact_genome::interface::alphabet::{Alphabet, AlphabetCharacter};

use super::MultialignMetric;

/// A pairwise metric with a pairwise scoring table.
pub struct PairwiseCostMetric<AlphabetType> {
    scoring_table: PairwiseScoringTable<AlphabetType>,
    character_counts: Vec<u8>,
    non_zero_character_counts: Vec<usize>,
    sequence_amount: i16,
    phantom_data: PhantomData<AlphabetType>,
}

struct PairwiseScoringTable<AlphabetType> {
    phantom_data: PhantomData<AlphabetType>,
}

impl<AlphabetType: Alphabet> PairwiseCostMetric<AlphabetType> {
    pub fn from_csv_file(path: impl AsRef<Path>, sequence_amount: usize) -> Result<Self> {
        Ok(Self {
            scoring_table: PairwiseScoringTable::from_csv_file(path)?,
            character_counts: vec![0; usize::from(AlphabetType::SIZE) + 1],
            non_zero_character_counts: Default::default(),
            // We multiply the i16 by itself later, so we restrict to i8 to make sure it does not overflow.
            sequence_amount: i8::try_from(sequence_amount)
                .with_context(|| format!("Metric supports at most {} sequences", i8::MAX))?
                .into(),
            phantom_data: PhantomData,
        })
    }
}

impl<AlphabetType: Alphabet> MultialignMetric<AlphabetType> for PairwiseCostMetric<AlphabetType> {
    fn reset_character_counts(&mut self) {
        self.character_counts.fill(0);
    }

    fn count_character(&mut self, character: &<AlphabetType as Alphabet>::CharacterType) {
        self.character_counts[usize::from(character.index())] += 1;
    }

    fn count_gap(&mut self) {
        self.character_counts[usize::from(AlphabetType::SIZE)] += 1;
    }

    fn compute_cost_increment<Cost: generic_a_star::cost::AStarCost>(&self) -> Cost
    where
        Cost::CostType: From<i16>,
    {
        todo!()
    }
}

impl<AlphabetType: Alphabet> PairwiseScoringTable<AlphabetType> {
    fn from_csv_file(path: impl AsRef<Path>) -> Result<Self> {
        todo!()
    }
}
