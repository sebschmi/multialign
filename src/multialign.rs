use std::{
    fmt::{Debug, Display},
    hash::Hash,
    marker::PhantomData,
};

use anyhow::Result;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::{cost::Cost, reset::Reset, AStar, AStarContext, AStarNode};
use log::info;

mod display;

trait NodeIdentifier: Debug + Display + Clone + Eq + Ord + Hash {
    fn create_root(sequence_amount: usize) -> Self;

    fn offset(&self, index: usize) -> usize;
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

    fn generate_successors(&mut self, _node: &Self::Node, _output: &mut impl Extend<Self::Node>) {
        todo!()
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
        64 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<64>>(
            sequences,
        ),
        65 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<65>>(
            sequences,
        ),
        66 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<66>>(
            sequences,
        ),
        67 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<67>>(
            sequences,
        ),
        68 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<68>>(
            sequences,
        ),
        69 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<69>>(
            sequences,
        ),
        70 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<70>>(
            sequences,
        ),
        71 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<71>>(
            sequences,
        ),
        72 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<72>>(
            sequences,
        ),
        73 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<73>>(
            sequences,
        ),
        74 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<74>>(
            sequences,
        ),
        75 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<75>>(
            sequences,
        ),
        76 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<76>>(
            sequences,
        ),
        77 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<77>>(
            sequences,
        ),
        78 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<78>>(
            sequences,
        ),
        79 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<79>>(
            sequences,
        ),
        80 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<80>>(
            sequences,
        ),
        81 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<81>>(
            sequences,
        ),
        82 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<82>>(
            sequences,
        ),
        83 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<83>>(
            sequences,
        ),
        84 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<84>>(
            sequences,
        ),
        85 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<85>>(
            sequences,
        ),
        86 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<86>>(
            sequences,
        ),
        87 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<87>>(
            sequences,
        ),
        88 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<88>>(
            sequences,
        ),
        89 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<89>>(
            sequences,
        ),
        90 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<90>>(
            sequences,
        ),
        91 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<91>>(
            sequences,
        ),
        92 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<92>>(
            sequences,
        ),
        93 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<93>>(
            sequences,
        ),
        94 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<94>>(
            sequences,
        ),
        95 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<95>>(
            sequences,
        ),
        96 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<96>>(
            sequences,
        ),
        97 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<97>>(
            sequences,
        ),
        98 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<98>>(
            sequences,
        ),
        99 => multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<99>>(
            sequences,
        ),
        100 => {
            multialign_astar_with_identifier::<AlphabetType, SequenceType, ArrayIdentifier<100>>(
                sequences,
            )
        }
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
    let mut a_star = AStar::new(Context::<_, _, Identifier>::new(sequences));
    match a_star.search() {
        generic_a_star::AStarResult::FoundTarget { cost, .. } => info!("Alignment cost {}", cost),
        generic_a_star::AStarResult::ExceededCostLimit { .. } => unreachable!("No cost limit set"),
        generic_a_star::AStarResult::ExceededMemoryLimit { .. } => {
            unreachable!("No memory limit set")
        }
        generic_a_star::AStarResult::NoTarget => unreachable!("Search always finds a target"),
    }

    Ok(())
}
