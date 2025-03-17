use std::{
    collections::HashSet,
    fmt::{Debug, Display},
    hash::Hash,
    marker::PhantomData,
    time::Instant,
    vec,
};

use anyhow::{bail, Context as _, Result};
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::{
    cost::{AStarCost, I16Cost},
    reset::Reset,
    AStar, AStarContext, AStarNode, AStarResult,
};
use log::info;
use metric::MultialignMetric;

mod display;
pub mod metric;

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
struct Node<Identifier: NodeIdentifier, Cost> {
    cost: Cost,
    identifier: Identifier,
    predecessor: Option<Identifier>,
}

impl<Identifier: NodeIdentifier, Cost: AStarCost> AStarNode for Node<Identifier, Cost> {
    type Identifier = Identifier;

    type EdgeType = Self;

    type Cost = Cost;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Self::Cost {
        self.cost
    }

    fn a_star_lower_bound(&self) -> Self::Cost {
        Self::Cost::zero()
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
    Cost,
    SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
    Identifier,
    Metric: MultialignMetric<AlphabetType>,
> {
    sequences: &'sequences [&'sequences SequenceType],
    metric: Metric,

    phantom_data: PhantomData<(Identifier, AlphabetType, Cost)>,
}

impl<
        AlphabetType: Alphabet,
        Cost: AStarCost,
        SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
        Identifier: NodeIdentifier,
        Metric: MultialignMetric<AlphabetType>,
    > AStarContext for Context<'_, AlphabetType, Cost, SequenceType, Identifier, Metric>
where
    Cost::CostType: From<i16>,
{
    type Node = Node<Identifier, Cost>;

    fn create_root(&self) -> Self::Node {
        Self::Node {
            cost: Cost::zero(),
            identifier: Identifier::create_root(self.sequences.len()),
            predecessor: None,
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        debug_assert!(self.sequences.len() < usize::BITS.try_into().unwrap());

        for gaps in 1..(1usize << self.sequences.len()) {
            // Compute next identifier and cost.
            let mut identifier = node.identifier.clone();
            self.metric.reset_character_counts();

            for (index, sequence) in self.sequences.iter().enumerate() {
                if gaps & (1 << index) != 0 && identifier.offset(index) < sequence.len() {
                    self.metric
                        .count_character(&sequence[identifier.offset(index)]);
                    identifier.increment(index);
                } else {
                    // Last entry represents a gap.
                    self.metric.count_gap();
                }
            }
            let identifier = identifier;

            // Compute cost increment.
            let cost_increment = self.metric.compute_cost_increment();

            output.extend(Some(Self::Node {
                cost: node.cost.checked_add(&cost_increment).unwrap(),
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
        Cost: AStarCost,
        SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
        Identifier: NodeIdentifier,
        Metric: MultialignMetric<AlphabetType>,
    > Reset for Context<'_, AlphabetType, Cost, SequenceType, Identifier, Metric>
{
    fn reset(&mut self) {
        // Do nothing.
    }
}

impl<
        'sequences,
        AlphabetType: Alphabet,
        Cost: AStarCost,
        SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
        Identifier: NodeIdentifier,
        Metric: MultialignMetric<AlphabetType>,
    > Context<'sequences, AlphabetType, Cost, SequenceType, Identifier, Metric>
{
    fn new(sequences: &'sequences [&'sequences SequenceType], metric: Metric) -> Self {
        Self {
            sequences,
            metric,
            phantom_data: PhantomData,
        }
    }
}

pub fn multialign_astar<
    AlphabetType: Alphabet + Debug + Clone + Eq + 'static,
    SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
    Metric: MultialignMetric<AlphabetType>,
>(
    sequences: &[&SequenceType],
    metric: Metric,
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
        2 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<2>, _>(
            sequences, metric,
        ),
        3 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<3>, _>(
            sequences, metric,
        ),
        4 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<4>, _>(
            sequences, metric,
        ),
        5 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<5>, _>(
            sequences, metric,
        ),
        6 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<6>, _>(
            sequences, metric,
        ),
        7 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<7>, _>(
            sequences, metric,
        ),
        8 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<8>, _>(
            sequences, metric,
        ),
        9 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<9>, _>(
            sequences, metric,
        ),
        10 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<10>, _>(
                sequences, metric,
            )
        }
        11 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<11>, _>(
                sequences, metric,
            )
        }
        12 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<12>, _>(
                sequences, metric,
            )
        }
        13 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<13>, _>(
                sequences, metric,
            )
        }
        14 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<14>, _>(
                sequences, metric,
            )
        }
        15 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<15>, _>(
                sequences, metric,
            )
        }
        16 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<16>, _>(
                sequences, metric,
            )
        }
        17 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<17>, _>(
                sequences, metric,
            )
        }
        18 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<18>, _>(
                sequences, metric,
            )
        }
        19 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<19>, _>(
                sequences, metric,
            )
        }
        20 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<20>, _>(
                sequences, metric,
            )
        }
        21 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<21>, _>(
                sequences, metric,
            )
        }
        22 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<22>, _>(
                sequences, metric,
            )
        }
        23 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<23>, _>(
                sequences, metric,
            )
        }
        24 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<24>, _>(
                sequences, metric,
            )
        }
        25 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<25>, _>(
                sequences, metric,
            )
        }
        26 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<26>, _>(
                sequences, metric,
            )
        }
        27 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<27>, _>(
                sequences, metric,
            )
        }
        28 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<28>, _>(
                sequences, metric,
            )
        }
        29 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<29>, _>(
                sequences, metric,
            )
        }
        30 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<30>, _>(
                sequences, metric,
            )
        }
        31 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<31>, _>(
                sequences, metric,
            )
        }
        32 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<32>, _>(
                sequences, metric,
            )
        }
        33 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<33>, _>(
                sequences, metric,
            )
        }
        34 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<34>, _>(
                sequences, metric,
            )
        }
        35 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<35>, _>(
                sequences, metric,
            )
        }
        36 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<36>, _>(
                sequences, metric,
            )
        }
        37 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<37>, _>(
                sequences, metric,
            )
        }
        38 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<38>, _>(
                sequences, metric,
            )
        }
        39 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<39>, _>(
                sequences, metric,
            )
        }
        40 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<40>, _>(
                sequences, metric,
            )
        }
        41 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<41>, _>(
                sequences, metric,
            )
        }
        42 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<42>, _>(
                sequences, metric,
            )
        }
        43 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<43>, _>(
                sequences, metric,
            )
        }
        44 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<44>, _>(
                sequences, metric,
            )
        }
        45 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<45>, _>(
                sequences, metric,
            )
        }
        46 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<46>, _>(
                sequences, metric,
            )
        }
        47 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<47>, _>(
                sequences, metric,
            )
        }
        48 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<48>, _>(
                sequences, metric,
            )
        }
        49 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<49>, _>(
                sequences, metric,
            )
        }
        50 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<50>, _>(
                sequences, metric,
            )
        }
        51 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<51>, _>(
                sequences, metric,
            )
        }
        52 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<52>, _>(
                sequences, metric,
            )
        }
        53 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<53>, _>(
                sequences, metric,
            )
        }
        54 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<54>, _>(
                sequences, metric,
            )
        }
        55 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<55>, _>(
                sequences, metric,
            )
        }
        56 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<56>, _>(
                sequences, metric,
            )
        }
        57 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<57>, _>(
                sequences, metric,
            )
        }
        58 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<58>, _>(
                sequences, metric,
            )
        }
        59 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<59>, _>(
                sequences, metric,
            )
        }
        60 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<60>, _>(
                sequences, metric,
            )
        }
        61 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<61>, _>(
                sequences, metric,
            )
        }
        62 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<62>, _>(
                sequences, metric,
            )
        }
        63 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<63>, _>(
                sequences, metric,
            )
        }
        _ => multialign_astar_with_identifier::<AlphabetType, SequenceType, VecIdentifier, _>(
            sequences, metric,
        ),
    }
}

fn multialign_astar_with_identifier<
    AlphabetType: Alphabet + Debug + Clone + Eq + 'static,
    SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
    Identifier: NodeIdentifier,
    Metric: MultialignMetric<AlphabetType>,
>(
    sequences: &[&SequenceType],
    metric: Metric,
) -> Result<()> {
    let start_time = Instant::now();
    let mut a_star = AStar::new(Context::<_, I16Cost, _, Identifier, _>::new(
        sequences, metric,
    ));
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
    Cost: AStarCost,
    SequenceType: GenomeSequence<AlphabetType, SequenceType> + ?Sized,
    Identifier: NodeIdentifier,
>(
    sequences: &[&SequenceType],
    edges: impl IntoIterator<Item = Node<Identifier, Cost>>,
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
