use std::{collections::BTreeMap, marker::PhantomData, path::Path};

use anyhow::{anyhow, ensure, Context, Result};
use compact_genome::interface::alphabet::{Alphabet, AlphabetCharacter};
use csv::ReaderBuilder;
use generic_a_star::cost::AStarCost;
use log::{info, trace};

use super::MultialignMetric;

/// A pairwise metric with a pairwise scoring table.
pub struct PairwiseCostMetric<AlphabetType> {
    cost_table: PairwiseCostTable<AlphabetType>,
    character_counts: Vec<u8>,
    non_zero_character_counts: Vec<usize>,
    phantom_data: PhantomData<AlphabetType>,
}

struct PairwiseCostTable<AlphabetType> {
    table: Vec<Option<i32>>,
    phantom_data: PhantomData<AlphabetType>,
}

impl<AlphabetType: Alphabet> PairwiseCostMetric<AlphabetType> {
    pub fn from_csv_file(path: impl AsRef<Path>) -> Result<Self> {
        Ok(Self {
            cost_table: PairwiseCostTable::from_csv_file(path)?,
            character_counts: vec![0; usize::from(AlphabetType::SIZE) + 1],
            non_zero_character_counts: Default::default(),
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

    fn compute_cost_increment<Cost: AStarCost>(&mut self) -> Result<Cost>
    where
        Cost::CostType: From<i32>,
    {
        let mut cost = Cost::zero();

        self.non_zero_character_counts.clear();
        for (index, count) in self.character_counts.iter().copied().enumerate() {
            if count == 0 {
                continue;
            }

            self.non_zero_character_counts.push(index);

            for other_index in self.non_zero_character_counts.iter().copied() {
                let other_count = self.character_counts[other_index];
                let count = i32::from(count);
                let other_count = i32::from(other_count);

                let multiplicity = count.checked_mul(other_count).unwrap();
                let base_cost = self.cost_table.cost(
                    if index == usize::from(AlphabetType::SIZE) {
                        None
                    } else {
                        Some(
                            AlphabetType::CharacterType::from_index(index.try_into().unwrap())
                                .unwrap(),
                        )
                    },
                    if other_index == usize::from(AlphabetType::SIZE) {
                        None
                    } else {
                        Some(
                            AlphabetType::CharacterType::from_index(
                                other_index.try_into().unwrap(),
                            )
                            .unwrap(),
                        )
                    },
                )?;
                cost += Cost::from(Cost::CostType::from(
                    multiplicity.checked_mul(base_cost).unwrap(),
                ));
            }
        }

        Ok(cost)
    }
}

impl<AlphabetType: Alphabet> PairwiseCostTable<AlphabetType> {
    fn from_csv_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        info!("Reading CSV file {path:?}");

        let mut reader = ReaderBuilder::new()
            .has_headers(false)
            .from_path(path)
            .with_context(|| format!("Error opening CSV file {path:?}"))?;
        let mut lines = reader.records().enumerate();
        let mut cost_map = BTreeMap::new();

        // Parse first line
        let first_line = lines
            .next()
            .map(|(_, first_line)| first_line)
            .ok_or_else(|| anyhow!("CSV file contains no lines"))?
            .with_context(|| "Error reading first CSV line")?;
        let mut character_to_column = vec![None; usize::from(AlphabetType::SIZE) + 1];
        let mut column_to_character = Vec::new();
        let gap_character_index = AlphabetType::SIZE;

        for (column, character) in first_line.iter().enumerate() {
            if column == 0 {
                ensure!(
                    character.trim().is_empty(),
                    "First column of first row must be empty, but was: {:?}",
                    character.trim()
                );
                column_to_character.push(None);
                continue;
            }

            let character = character.trim();
            ensure!(character.chars().count() == 1, "First row must contain a single character in each column, except for the first column which must be empty");
            let character = character.chars().next().unwrap();

            if character == '*' || character == '-' {
                ensure!(
                    character_to_column[usize::from(gap_character_index)].is_none(),
                    "First row contained a character twice: {character}"
                );
                character_to_column[usize::from(gap_character_index)] = Some(column);
                column_to_character.push(Some(gap_character_index));
                continue;
            }

            let character = character.try_into().with_context(|| {
                "Character must be a valid ASCII character, but is {character:?}"
            })?;
            let character = AlphabetType::ascii_to_character(character)
                .with_context(|| "Character must be a valid alphabet character")?;
            let index = usize::from(character.index());

            ensure!(
                character_to_column[index].is_none(),
                "First row contained a character twice: {character}"
            );
            character_to_column[index] = Some(column);
            column_to_character.push(Some(character.index()));
        }

        // Parse further lines
        let mut character_to_row = vec![None; usize::from(AlphabetType::SIZE) + 1];
        let mut row_to_character = vec![None];
        for (row, line) in lines {
            let line = line.with_context(|| "Error reading CSV line")?;

            for (column, cost) in line.iter().enumerate() {
                if column == 0 {
                    let character = cost.trim();
                    ensure!(
                        character.chars().count() == 1,
                        "First column of rows after the first must contain a single character"
                    );
                    let character = character.chars().next().unwrap();

                    if character == '*' || character == '-' {
                        ensure!(
                            character_to_row[usize::from(gap_character_index)].is_none(),
                            "First column contained a character twice: {character}"
                        );
                        character_to_row[usize::from(gap_character_index)] = Some(row);
                        row_to_character.push(Some(gap_character_index));
                        continue;
                    }

                    let character = character.try_into().with_context(|| {
                        "Character must be a valid ASCII character, but is {character:?}"
                    })?;
                    let character = AlphabetType::ascii_to_character(character)
                        .with_context(|| "Character must be a valid alphabet character")?;
                    let index = usize::from(character.index());

                    ensure!(
                        character_to_row[index].is_none(),
                        "First column contained a character twice: {character}"
                    );
                    character_to_row[index] = Some(column);
                    row_to_character.push(Some(character.index()));
                    continue;
                }

                let cost = cost.trim();
                let cost = cost
                    .parse()
                    .with_context(|| format!("Error parsing '{cost}' as i32"))?;

                let from = row_to_character[row].unwrap();
                let to = column_to_character[column].ok_or_else(|| {
                    anyhow!("Subsequent row contains more columns than the first row")
                })?;
                let from = if from == gap_character_index {
                    None
                } else {
                    Some(AlphabetType::CharacterType::from_index(from).unwrap())
                };
                let to = if to == gap_character_index {
                    None
                } else {
                    Some(AlphabetType::CharacterType::from_index(to).unwrap())
                };
                let previous_cost = cost_map.insert((from, to), cost);
                debug_assert!(previous_cost.is_none());
            }
        }

        // Transform map into table
        let mut table = Vec::with_capacity((usize::from(AlphabetType::SIZE) + 1) << 1);
        for from in AlphabetType::iter().map(Some).chain([None]) {
            for to in AlphabetType::iter().map(Some).chain([None]) {
                trace!(
                    "from_index: {:?}; to_index: {:?}; gap_character_index: {}",
                    from.as_ref().map(AlphabetCharacter::index),
                    to.as_ref().map(AlphabetCharacter::index),
                    gap_character_index,
                );
                let cost = cost_map.get(&(from.clone(), to.clone())).copied();
                table.push(cost);

                let from_index = from
                    .as_ref()
                    .map(AlphabetCharacter::index)
                    .unwrap_or(gap_character_index);
                let to_index = to
                    .as_ref()
                    .map(AlphabetCharacter::index)
                    .unwrap_or(gap_character_index);

                // Ensure symmetry
                if to_index < from_index {
                    let other_cost = table[usize::from(to_index)
                        * (usize::from(AlphabetType::SIZE) + 1)
                        + usize::from(from_index)];
                    ensure!(
                        other_cost == cost,
                        "Assymetric entry found for row {}, column {}: {:?} != {:?}",
                        if let Some(from) = from {
                            from.to_string()
                        } else {
                            "gap".to_string()
                        },
                        if let Some(to) = to {
                            to.to_string()
                        } else {
                            "gap".to_string()
                        },
                        cost,
                        other_cost,
                    )
                }
            }
        }

        Ok(Self {
            table,
            phantom_data: PhantomData,
        })
    }

    fn cost(
        &self,
        from: Option<AlphabetType::CharacterType>,
        to: Option<AlphabetType::CharacterType>,
    ) -> Result<i32> {
        let gap_character_index = usize::from(AlphabetType::SIZE);
        let from = usize::from(
            from.as_ref()
                .map(AlphabetCharacter::index)
                .unwrap_or(AlphabetType::SIZE),
        );
        let to = usize::from(
            to.as_ref()
                .map(AlphabetCharacter::index)
                .unwrap_or(AlphabetType::SIZE),
        );

        self.table[from * (usize::from(AlphabetType::SIZE) + 1) + to].ok_or_else(|| {
            anyhow!(
                "Missing cost for row {}, column {}",
                if from != gap_character_index {
                    from.to_string()
                } else {
                    "gap".to_string()
                },
                if to != gap_character_index {
                    to.to_string()
                } else {
                    "gap".to_string()
                }
            )
        })
    }
}
