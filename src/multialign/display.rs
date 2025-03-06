use std::fmt::Display;

use super::{ArrayIdentifier, Node, NodeIdentifier, VecIdentifier};

impl<const SEQUENCE_AMOUNT: usize> Display for ArrayIdentifier<SEQUENCE_AMOUNT> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "(")?;

        let mut once = true;
        for offset in self.offsets {
            if once {
                once = false;
            } else {
                write!(f, ", ")?;
            }

            write!(f, "{offset}")?;
        }

        write!(f, ")")
    }
}

impl Display for VecIdentifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "(")?;

        let mut once = true;
        for offset in &self.offsets {
            if once {
                once = false;
            } else {
                write!(f, ", ")?;
            }

            write!(f, "{offset}")?;
        }

        write!(f, ")")
    }
}

impl<Identifier: NodeIdentifier, Cost: Display> Display for Node<Identifier, Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Node[{} -> {}; {}]",
            self.predecessor
                .as_ref()
                .map(ToString::to_string)
                .unwrap_or_else(|| "None".to_string()),
            self.identifier,
            self.cost
        )
    }
}
