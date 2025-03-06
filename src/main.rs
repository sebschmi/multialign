use std::{fmt::Debug, path::PathBuf};

use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use compact_genome::{
    implementation::{
        alphabets::{
            dna_alphabet::DnaAlphabet, dna_alphabet_or_n::DnaAlphabetOrN,
            dna_iupac_nucleic_acid_alphabet::DnaIupacNucleicAcidAlphabet,
            famsa_amino_acid_alphabet::FamsaAminoAcidAlphabet,
            iupac_amino_acid_alphabet::IupacAminoAcidAlphabet, rna_alphabet::RnaAlphabet,
            rna_alphabet_or_n::RnaAlphabetOrN,
            rna_iupac_nucleic_acid_alphabet::RnaIupacNucleicAcidAlphabet,
        },
        DefaultSequenceStore,
    },
    interface::{alphabet::Alphabet, sequence::GenomeSequence, sequence_store::SequenceStore},
    io::fasta::read_fasta_file,
};
use log::{error, info, LevelFilter};
use multialign::multialign_astar;
use simplelog::{ColorChoice, TermLogger, TerminalMode};

mod multialign;

#[derive(Parser)]
struct Cli {
    /// The minimum importance of log messages to output.
    #[clap(long, short = 'l', default_value = "info")]
    log_level: LevelFilter,

    /// The input sequences.
    /// This can be multiple fasta files where each file may contain multiple sequences as fasta records.
    #[clap(long, short = 'i')]
    input: Vec<PathBuf>,

    /// The alphabet present in the input files.
    ///
    /// This must also match the alphabet used in the config.
    #[clap(long, short = 'a', default_value = "famsa-amino-acid")]
    alphabet: InputAlphabet,

    /// A string of (ASCII) characters that should be skipped in the input fasta.
    ///
    /// For example, `-` characters caused by alignment hints can be skipped this way.
    #[clap(long, default_value = "")]
    skip_characters: String,
}

#[derive(Debug, Clone, Eq, PartialEq, ValueEnum)]
enum InputAlphabet {
    Dna,
    DnaN,
    Rna,
    RnaN,
    DnaIupac,
    RnaIupac,
    /// The IUPAC amino acid alphabet.
    IupacAminoAcid,
    /// The FAMSA amino acid alphabet.
    FamsaAminoAcid,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    TermLogger::init(
        cli.log_level,
        Default::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )?;

    info!("Logging initialised");

    match cli.alphabet {
        InputAlphabet::Dna => execute_with_alphabet::<DnaAlphabet>(cli),
        InputAlphabet::DnaN => execute_with_alphabet::<DnaAlphabetOrN>(cli),
        InputAlphabet::Rna => execute_with_alphabet::<RnaAlphabet>(cli),
        InputAlphabet::RnaN => execute_with_alphabet::<RnaAlphabetOrN>(cli),
        InputAlphabet::DnaIupac => execute_with_alphabet::<DnaIupacNucleicAcidAlphabet>(cli),
        InputAlphabet::RnaIupac => execute_with_alphabet::<RnaIupacNucleicAcidAlphabet>(cli),
        InputAlphabet::IupacAminoAcid => execute_with_alphabet::<IupacAminoAcidAlphabet>(cli),
        InputAlphabet::FamsaAminoAcid => execute_with_alphabet::<FamsaAminoAcidAlphabet>(cli),
    }?;

    info!("Terminating");

    Ok(())
}

fn execute_with_alphabet<AlphabetType: Alphabet + Debug + Clone + Eq + 'static>(
    cli: Cli,
) -> Result<()> {
    if cli.input.is_empty() {
        bail!("No input files given");
    }

    let mut skip_characters = Vec::new();
    for character in cli.skip_characters.bytes().map(usize::from) {
        if skip_characters.len() <= character {
            skip_characters.resize(character + 1, false);
        }
        skip_characters[character] = true;
    }
    let skip_characters = skip_characters;

    let mut sequence_store = DefaultSequenceStore::<AlphabetType>::new();
    let mut records = Vec::new();
    for path in &cli.input {
        info!("Loading fasta file {path:?}");
        let path_records =
            read_fasta_file(path, &mut sequence_store, false, true, &skip_characters)
                .with_context(|| format!("Error loading file: {path:?}"))?;

        for mut record in path_records {
            if cli.input.len() > 1 {
                record.id = format!("{path:?}-{}", record.id);
            }

            records.push(record);
        }
    }

    if records.is_empty() {
        bail!("Found no fasta records in input files");
    } else if records.len() == 1 {
        bail!("Found only one fasta record in input files");
    }

    let mut record_ids: Vec<_> = records.iter().map(|record| record.id.clone()).collect();
    record_ids.sort_unstable();
    let duplicate_ids = list_duplicates(&record_ids);
    if !duplicate_ids.is_empty() {
        for duplicate_id in &duplicate_ids {
            error!("Found duplicate id {duplicate_id}");
        }

        bail!("Found {} distinct duplicate ids", duplicate_ids.len());
    }

    info!("Loaded {} sequences", records.len());

    let sequences: Vec<_> = records
        .iter()
        .map(|record| {
            sequence_store
                .get(&record.sequence_handle)
                .as_genome_subsequence()
        })
        .collect();
    multialign_astar(&sequences)
}

fn list_duplicates<T: Eq + Ord>(slice: &[T]) -> Vec<&T> {
    debug_assert!(slice.is_sorted());

    if slice.is_empty() {
        return Default::default();
    }

    let mut duplicates = Vec::new();
    let mut last_element = &slice[0];

    for element in &slice[1..] {
        if element == last_element && duplicates.last() != Some(&element) {
            duplicates.push(element);
        }

        last_element = element;
    }

    duplicates
}
