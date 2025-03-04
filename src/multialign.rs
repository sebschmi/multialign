use std::fmt::Debug;

use anyhow::Result;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};

pub fn multialign_astar<
    AlphabetType: Alphabet + Debug + Clone + Eq + 'static,
    Genome: GenomeSequence<AlphabetType, Genome> + ?Sized,
>(
    _sequences: &[&Genome],
) -> Result<()> {
    Ok(())
}
