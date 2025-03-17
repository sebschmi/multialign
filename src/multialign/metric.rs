use anyhow::Result;
use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::AStarCost;

pub mod pairwise_cost_metric;
pub mod pairwise_match_metric;

pub trait MultialignMetric<AlphabetType: Alphabet> {
    fn reset_character_counts(&mut self);

    fn count_character(&mut self, character: &AlphabetType::CharacterType);

    fn count_gap(&mut self);

    fn compute_cost_increment<Cost: AStarCost>(&mut self) -> Result<Cost>
    where
        Cost::CostType: From<i32>;
}
