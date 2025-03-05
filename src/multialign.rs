use std::{
    collections::HashSet,
    fmt::{Debug, Display},
    hash::Hash,
    marker::PhantomData,
    time::Instant,
    vec,
};

use anyhow::{bail, Context as _, Result};
use compact_genome::interface::{
    alphabet::{Alphabet, AlphabetCharacter},
    sequence::GenomeSequence,
};
use generic_a_star::{cost::Cost, reset::Reset, AStar, AStarContext, AStarNode, AStarResult};
use log::info;

mod display;

trait NodeIdentifier: Debug + Display + Clone + Eq + Ord + Hash {
    fn create_root(sequence_amount: usize) -> Self;

    fn offset(&self, index: usize) -> usize;

    fn increment(&mut self, index: usize);
}

#[derive(Debug, Clone, Hash, Eq, PartialEq, Ord, PartialOrd)]
struct ArrayIdentifier<const SEQUENCE_AMOUNT: usize> {
    offsets: [usize; SEQUENCE_AMOUNT],
}

#[derive(Debug, Clone, Hash, Eq, PartialEq, Ord, PartialOrd)]
struct VecIdentifier {
    offsets: Vec<usize>,
}

impl<const SEQUENCE_AMOUNT: usize> NodeIdentifier for ArrayIdentifier<SEQUENCE_AMOUNT> {
    fn create_root(sequence_amount: usize) -> Self {
        assert_eq!(sequence_amount, SEQUENCE_AMOUNT);
        Self {
            offsets: [0; SEQUENCE_AMOUNT],
        }
    }

    fn offset(&self, index: usize) -> usize {
        self.offsets[index]
    }

    fn increment(&mut self, index: usize) {
        self.offsets[index] += 1;
    }
}

impl NodeIdentifier for VecIdentifier {
    fn create_root(sequence_amount: usize) -> Self {
        Self {
            offsets: vec![0; sequence_amount],
        }
    }

    fn offset(&self, index: usize) -> usize {
        self.offsets[index]
    }

    fn increment(&mut self, index: usize) {
        self.offsets[index] += 1;
    }
}

#[derive(Debug, Clone, Hash, Eq, PartialEq, Ord, PartialOrd)]
struct Node<Identifier: NodeIdentifier> {
    cost: Cost,
    identifier: Identifier,
    predecessor: Option<Identifier>,
}

impl<Identifier: NodeIdentifier> AStarNode for Node<Identifier> {
    type Identifier = Identifier;

    type EdgeType = Self;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Cost {
        self.cost
    }

    fn a_star_lower_bound(&self) -> Cost {
        Cost::ZERO
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        self.predecessor.as_ref()
    }

    fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
        self.predecessor.as_ref().map(|_| self.clone())
    }
}

struct Context<
    'sequences,
    AlphabetType: Alphabet,
    SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
    Identifier: NodeIdentifier,
> {
    sequences: &'sequences [&'sequences SequenceType],
    character_counts: Vec<u8>,

    phantom_data: PhantomData<(Identifier, AlphabetType)>,
}

impl<
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
        Identifier: NodeIdentifier,
    > AStarContext for Context<'_, AlphabetType, SequenceType, Identifier>
{
    type Node = Node<Identifier>;

    fn create_root(&self) -> Self::Node {
        Self::Node {
            cost: Cost::ZERO,
            identifier: Identifier::create_root(self.sequences.len()),
            predecessor: None,
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        for gaps in 1..(1usize << self.sequences.len()) {
            // Compute next identifier and cost.
            let mut identifier = node.identifier.clone();
            self.character_counts.fill(0);

            for (index, sequence) in self.sequences.iter().enumerate() {
                if gaps & (1 << index) != 0 && identifier.offset(index) < sequence.len() {
                    self.character_counts
                        [usize::from(sequence[identifier.offset(index)].index())] += 1;
                    identifier.increment(index);
                } else {
                    // Last entry represents a gap.
                    self.character_counts[usize::from(AlphabetType::SIZE)] += 1;
                }
            }
            let identifier = identifier;

            // Compute cost increment.
            let score_increment =
                self.character_counts
                    .iter()
                    .fold(Cost::ZERO, |score, character_count| {
                        let character_count = u64::from(*character_count);
                        let character_score = if character_count >= 2 {
                            Cost::from(character_count * (character_count - 1) / 2)
                        } else {
                            Cost::ZERO
                        };
                        score + character_score
                    });
            let max_score =
                Cost::from((self.sequences.len() * (self.sequences.len() - 1) / 2) as u64);
            let cost_increment = max_score - score_increment;

            output.extend(Some(Self::Node {
                cost: node.cost + cost_increment,
                identifier,
                predecessor: Some(node.identifier.clone()),
            }));
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        self.sequences.iter().enumerate().all(|(index, sequence)| {
            debug_assert!(node.identifier.offset(index) <= sequence.len());
            node.identifier.offset(index) == sequence.len()
        })
    }

    fn cost_limit(&self) -> Option<Cost> {
        None
    }

    fn memory_limit(&self) -> Option<usize> {
        None
    }
}

impl<
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
        Identifier: NodeIdentifier,
    > Reset for Context<'_, AlphabetType, SequenceType, Identifier>
{
    fn reset(&mut self) {
        // Do nothing.
    }
}

impl<
        'sequences,
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
        Identifier: NodeIdentifier,
    > Context<'sequences, AlphabetType, SequenceType, Identifier>
{
    fn new(sequences: &'sequences [&'sequences SequenceType]) -> Self {
        Self {
            sequences,
            character_counts: vec![0; usize::from(AlphabetType::SIZE) + 1],
            phantom_data: PhantomData,
        }
    }
}

pub fn multialign_astar<
    AlphabetType: Alphabet + Debug + Clone + Eq + 'static,
    SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
>(
    sequences: &[&SequenceType],
) -> Result<()> {
    info!("Aligning {} sequences", sequences.len());

    let max_sequence_amount = usize::BITS - 1;
    let sequence_len_u32: u32 = sequences.len().try_into().with_context(|| {
        format!(
            "Exceeded maximum supported sequence amount: {} > {}",
            sequences.len(),
            max_sequence_amount
        )
    })?;
    if sequence_len_u32 > max_sequence_amount {
        bail!(
            "Exceeded maximum supported sequence amount: {} > {}",
            sequences.len(),
            max_sequence_amount
        );
    }

    match sequences.len() {
        0 | 1 => panic!("Called multialign_astar with less than two sequences"),
        2 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<2>>(
            sequences,
        ),
        3 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<3>>(
            sequences,
        ),
        4 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<4>>(
            sequences,
        ),
        5 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<5>>(
            sequences,
        ),
        6 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<6>>(
            sequences,
        ),
        7 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<7>>(
            sequences,
        ),
        8 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<8>>(
            sequences,
        ),
        9 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<9>>(
            sequences,
        ),
        10 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<10>>(
            sequences,
        ),
        11 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<11>>(
            sequences,
        ),
        12 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<12>>(
            sequences,
        ),
        13 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<13>>(
            sequences,
        ),
        14 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<14>>(
            sequences,
        ),
        15 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<15>>(
            sequences,
        ),
        16 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<16>>(
            sequences,
        ),
        17 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<17>>(
            sequences,
        ),
        18 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<18>>(
            sequences,
        ),
        19 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<19>>(
            sequences,
        ),
        20 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<20>>(
            sequences,
        ),
        21 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<21>>(
            sequences,
        ),
        22 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<22>>(
            sequences,
        ),
        23 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<23>>(
            sequences,
        ),
        24 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<24>>(
            sequences,
        ),
        25 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<25>>(
            sequences,
        ),
        26 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<26>>(
            sequences,
        ),
        27 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<27>>(
            sequences,
        ),
        28 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<28>>(
            sequences,
        ),
        29 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<29>>(
            sequences,
        ),
        30 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<30>>(
            sequences,
        ),
        31 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<31>>(
            sequences,
        ),
        32 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<32>>(
            sequences,
        ),
        33 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<33>>(
            sequences,
        ),
        34 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<34>>(
            sequences,
        ),
        35 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<35>>(
            sequences,
        ),
        36 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<36>>(
            sequences,
        ),
        37 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<37>>(
            sequences,
        ),
        38 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<38>>(
            sequences,
        ),
        39 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<39>>(
            sequences,
        ),
        40 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<40>>(
            sequences,
        ),
        41 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<41>>(
            sequences,
        ),
        42 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<42>>(
            sequences,
        ),
        43 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<43>>(
            sequences,
        ),
        44 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<44>>(
            sequences,
        ),
        45 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<45>>(
            sequences,
        ),
        46 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<46>>(
            sequences,
        ),
        47 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<47>>(
            sequences,
        ),
        48 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<48>>(
            sequences,
        ),
        49 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<49>>(
            sequences,
        ),
        50 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<50>>(
            sequences,
        ),
        51 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<51>>(
            sequences,
        ),
        52 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<52>>(
            sequences,
        ),
        53 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<53>>(
            sequences,
        ),
        54 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<54>>(
            sequences,
        ),
        55 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<55>>(
            sequences,
        ),
        56 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<56>>(
            sequences,
        ),
        57 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<57>>(
            sequences,
        ),
        58 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<58>>(
            sequences,
        ),
        59 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<59>>(
            sequences,
        ),
        60 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<60>>(
            sequences,
        ),
        61 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<61>>(
            sequences,
        ),
        62 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<62>>(
            sequences,
        ),
        63 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<63>>(
            sequences,
        ),
        _ => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, VecIdentifier>(sequences)
        }
    }
}

fn multialign_astar_with_identifier<
    AlphabetType: Alphabet + Debug + Clone + Eq + 'static,
    SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
    Identifier: NodeIdentifier,
>(
    sequences: &[&SequenceType],
) -> Result<()> {
    let start_time = Instant::now();
    let mut a_star = AStar::new(Context::<_, _, Identifier>::new(sequences));
    a_star.initialise();

    match a_star.search() {
        AStarResult::FoundTarget { cost, .. } => info!("Alignment cost {}", cost),
        AStarResult::ExceededCostLimit { .. } => unreachable!("No cost limit set"),
        AStarResult::ExceededMemoryLimit { .. } => {
            unreachable!("No memory limit set")
        }
        AStarResult::NoTarget => unreachable!("Search always finds a target"),
    }

    let end_time = Instant::now();
    let duration = end_time - start_time;

    info!("Runtime: {:.2}s", duration.as_secs_f64());
    info!("Performance: {:?}", a_star.performance_counters());
    info!(
        "Alignment: {}",
        backtrack_cigar(sequences, a_star.backtrack())
    );

    Ok(())
}

fn backtrack_cigar<
    AlphabetType: Alphabet + Debug + Clone + Eq + 'static,
    SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
    Identifier: NodeIdentifier,
>(
    sequences: &[&SequenceType],
    edges: impl IntoIterator<Item = Node<Identifier>>,
) -> String {
    enum CigarElement {
        Match { amount: usize },
        Mismatch { column: Vec<Option<char>> },
    }

    let mut cigar = Vec::new();

    for edge in edges {
        let mut column = Vec::new();

        for (index, sequence) in sequences.iter().enumerate() {
            let predecessor_offset = edge.predecessor.as_ref().unwrap().offset(index);
            let offset = edge.identifier.offset(index);

            if predecessor_offset == offset {
                column.push(None);
            } else {
                debug_assert_eq!(predecessor_offset + 1, offset);
                column.push(Some(sequence[predecessor_offset].clone().into()));
            }
        }

        let column_set: HashSet<_> = column.iter().copied().collect();
        if column_set.len() == 1 {
            if let Some(CigarElement::Match { amount }) = cigar.last_mut() {
                *amount += 1;
            } else {
                cigar.push(CigarElement::Match { amount: 1 });
            }
        } else {
            cigar.push(CigarElement::Mismatch { column });
        }
    }

    let mut cigar_string = String::new();
    for element in cigar.iter().rev() {
        match element {
            CigarElement::Match { amount } => cigar_string.push_str(&format!("{amount}M")),
            CigarElement::Mismatch { column } => {
                cigar_string.push('[');
                for character in column {
                    cigar_string.push(character.unwrap_or('-'));
                }
                cigar_string.push(']');
            }
        }
    }

    cigar_string
}
