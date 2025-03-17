use std::marker::PhantomData;

use anyhow::{Context, Result};
use compact_genome::interface::alphabet::{Alphabet, AlphabetCharacter};

use super::MultialignMetric;

/// A pairwise metric that scores matches with zero and everything else with one.
///
/// Specifically, pairs of gaps are scored with zero as well.
pub struct PairwiseMatchMetric<AlphabetType: Alphabet> {
    character_counts: Vec<u8>,
    sequence_amount: i32,
    phantom_data: PhantomData<AlphabetType>,
}

impl<AlphabetType: Alphabet> PairwiseMatchMetric<AlphabetType> {
    pub fn new(sequence_amount: usize) -> Result<Self> {
        Ok(Self {
            character_counts: vec![0; usize::from(AlphabetType::SIZE) + 1],
            // We multiply the i32 by itself later, so we restrict to i8 to make sure it does not overflow.
            sequence_amount: i8::try_from(sequence_amount)
                .with_context(|| format!("Metric supports at most {} sequences", i8::MAX))?
                .into(),
            phantom_data: PhantomData,
        })
    }
}

impl<AlphabetType: Alphabet> MultialignMetric<AlphabetType> for PairwiseMatchMetric<AlphabetType> {
    fn reset_character_counts(&mut self) {
        self.character_counts.fill(0);
    }

    fn count_character(&mut self, character: &<AlphabetType as Alphabet>::CharacterType) {
        self.character_counts[usize::from(character.index())] += 1;
    }

    fn count_gap(&mut self) {
        self.character_counts[usize::from(AlphabetType::SIZE)] += 1;
    }

    fn compute_cost_increment<Cost: generic_a_star::cost::AStarCost>(&mut self) -> Result<Cost>
    where
        Cost::CostType: From<i32>,
    {
        let score_increment =
            self.character_counts
                .iter()
                .fold(Cost::zero(), |score, character_count| {
                    let character_count = i32::from(*character_count);
                    let character_score = if character_count >= 2 {
                        Cost::from(Cost::CostType::from(
                            character_count
                                .checked_mul(character_count.checked_sub(1).unwrap())
                                .unwrap()
                                .checked_div(2)
                                .unwrap(),
                        ))
                    } else {
                        Cost::zero()
                    };

                    score.checked_add(&character_score).unwrap()
                });
        let max_score = Cost::from(Cost::CostType::from(
            self.sequence_amount
                .checked_mul(self.sequence_amount.checked_sub(1).unwrap())
                .unwrap()
                .checked_div(2)
                .unwrap(),
        ));
        Ok(max_score.checked_sub(&score_increment).unwrap())
    }
}
