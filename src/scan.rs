use crate::kr_hash::KarpRabinHash;

use ahash::RandomState;
use seq_io::fasta;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

type BString = Vec<u8>;

const END_OF_DICT: u8 = 0;
const END_OF_WORD: u8 = 1;
const DOLLAR: u8 = 2;

struct WordStatistics {
    word: BString,
    occurence: u32,
    rank: u32,
}

impl Default for WordStatistics {
    fn default() -> Self {
        Self {
            word: vec![],
            occurence: 0,
            rank: 0,
        }
    }
}

struct TempFiles {
    files: Vec<BufWriter<File>>,
}

impl TempFiles {
    fn new(basename: &Path) -> TempFiles {
        let files = [".parse", ".last", ".dict", ".occ"]
            .iter()
            .map(|e| TempFiles::open_file(basename.with_extension(e)))
            .collect::<Vec<_>>();

        TempFiles { files }
    }

    fn write_to_parse_file(&mut self, data: &[u8]) {
        self.files[0].write(data).unwrap();
    }

    fn write_to_last_file(&mut self, data: &[u8]) {
        self.files[1].write(data).unwrap();
    }

    fn write_to_dict_file(&mut self, data: &[u8]) {
        self.files[2].write(data).unwrap();
    }

    fn write_to_occ_file(&mut self, data: &[u8]) {
        self.files[3].write(data).unwrap();
    }

    fn open_file(filename: impl AsRef<Path>) -> BufWriter<File> {
        File::create(filename)
            .map(BufWriter::new)
            .expect("Unable to create file")
    }
}

pub struct Scanner {
    map: HashMap<u64, WordStatistics, RandomState>,
    hashes: Vec<u64>,
    tmp_files: TempFiles,
}

impl Scanner {
    pub fn new(basename: &Path) -> Self {
        Scanner {
            map: HashMap::default(),
            hashes: Vec::new(),
            tmp_files: TempFiles::new(basename),
        }
    }

    fn save_or_update_word(&mut self, word: &[u8]) {
        let hash = KarpRabinHash::kr_hash(word);

        self.hashes.push(hash);

        let entry = self.map.entry(hash).or_default();

        if entry.occurence == 0 {
            entry.word = word.to_vec();
        } else if entry.word != word {
            eprintln!("Hash collision detected\n")
        }
        entry.occurence = entry
            .occurence
            .checked_add(1)
            .expect("Maximum occurence of dictionary word exceeded");
    }

    pub fn serialize(&mut self) {
        let mut values = self.map.values_mut().collect::<Vec<_>>();
        values.sort_by(|a, b| a.word.cmp(&b.word));

        let mut word_rank = 1;

        for v in values {
            self.tmp_files.write_to_dict_file(&v.word);
            self.tmp_files.write_to_dict_file(&[END_OF_WORD]);
            self.tmp_files.write_to_occ_file(&v.occurence.to_ne_bytes());
            (v.rank, word_rank) = (word_rank, word_rank + 1);
        }

        self.tmp_files.write_to_dict_file(&[END_OF_DICT]);
        for hash in &self.hashes {
            self.tmp_files
                .write_to_parse_file(&self.map[hash].rank.to_ne_bytes());
        }
    }

    pub fn process_fasta<T: Read>(
        &mut self,
        mut fasta: fasta::Reader<T>,
        window: u32,
        phrase_modulo: u32,
    ) {
        let mut word = BString::new();
        let mut kr_hash = KarpRabinHash::new(window);
        // let mut reader: Box<Reader<dyn Read>> = if filename.into().is_empty() {
        //     Box::new(fasta::Reader::from_path(filename).expect("Unable to open FASTA file"));
        // } else {
        //     let stdin = std::io::stdin();
        //     Box::new(fasta::Reader::new(stdin.lock()));
        // };

        word.push(DOLLAR);
        let mut process_line = |line: &[u8]| {
            for &c in line {
                word.push(c);
                let hash = kr_hash.add_character(c);
                if hash % phrase_modulo as u64 == 0 && word.len() > window as usize {
                    self.save_or_update_word(&word);
                    let overlap = word.len() - window as usize;
                    self.tmp_files
                        .write_to_last_file(&word[overlap - 1..overlap]);
                    word.drain(..overlap);
                }
            }
        };

        while let Some(result) = fasta.next() {
            let record = result.unwrap();

            record.seq_lines().for_each(&mut process_line);
        }

        for _ in 0..window {
            word.push(DOLLAR);
        }
        self.save_or_update_word(&word);
        let overlap = word.len() - window as usize;
        self.tmp_files
            .write_to_last_file(&word[overlap - 1..overlap]);
    }
}
