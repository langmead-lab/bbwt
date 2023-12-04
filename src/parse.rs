use crate::gsacak_c::{gsacak_u8, sacak_u32};

use byteorder::{NativeEndian, ReadBytesExt, WriteBytesExt};

use std::cmp::{Ordering, PartialEq, PartialOrd, Reverse};
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

const END_OF_WORD: u8 = 1;

fn read_file_into_vec_u8<P: AsRef<Path>>(path: P, vec: &mut Vec<u8>) {
    let mut f = File::open(path).expect("Unable to open file");
    let size = f.metadata().unwrap().len() as usize;

    if vec.len() == 0 {
        vec.resize(size + 1, 0);
    }

    f.read(&mut vec[0..size]).unwrap();
}

fn read_file_into_vec_u32<P: AsRef<Path>>(path: P, vec: &mut Vec<u32>) {
    let mut f = File::open(path).expect("Unable to open file");
    let num_u32s = (f.metadata().unwrap().len() / 4) as usize;

    if vec.len() == 0 {
        vec.resize(num_u32s + 1, 0);
    }

    f.read_u32_into::<NativeEndian>(&mut vec[0..num_u32s as usize])
        .unwrap();
}

pub struct Parse<'a> {
    parse: Vec<u32>,
    bwt: Vec<u32>,
    basename: &'a Path,
}

impl<'a> Parse<'a> {
    pub fn new(basename: &'a Path) -> Self {
        Parse {
            parse: vec![],
            bwt: vec![],
            basename,
        }
    }

    fn read_parse(&mut self) {
        read_file_into_vec_u32(self.basename.with_extension(".parse"), &mut self.parse);
    }

    fn compute_suffix_array(&mut self) {
        self.bwt.resize(self.parse.len(), 0);

        unsafe {
            sacak_u32(&self.parse, &mut self.bwt);
        }
    }

    pub fn transform_text_into_bwt(&mut self) {
        let mut output = File::create(self.basename.with_extension(".bwlast"))
            .map(BufWriter::new)
            .unwrap();

        self.read_parse();
        let n = self.parse.len() - 1;

        self.compute_suffix_array();
        {
            self.bwt[0] = self.parse[n - 1];
            let mut last = vec![0; n];
            read_file_into_vec_u8(self.basename.with_extension(".last"), &mut last);
            output.write(&[last[n - 2]]).unwrap();

            for i in 1..=n {
                if self.bwt[i] == 0 {
                    output.write(&[0]).unwrap();
                    continue;
                }
                if self.bwt[i] == 1 {
                    output.write(&[last[n - 1]]).unwrap();
                } else {
                    output.write(&[last[self.bwt[i] as usize - 2]]).unwrap();
                }
                self.bwt[i] = self.parse[self.bwt[i] as usize - 1];
            }
        }
    }

    pub fn contruct_inverse_list(&mut self) {
        let n = self.parse.len() - 1;
        let &alphabet_size = self.parse[0..n].iter().max().unwrap();
        let mut occ = vec![0; alphabet_size as usize + 1];

        read_file_into_vec_u32(self.basename.with_extension(".occ"), &mut occ);
        occ.rotate_right(1);
        occ[0] = 1;

        let mut frequency_vector = vec![0u32; alphabet_size as usize + 1];
        for i in 1..=alphabet_size as usize {
            frequency_vector[i] = frequency_vector[i - 1] + occ[i - 1];
        }

        // reuse the parse to compute store the inverse list
        let n = self.parse.len();
        let ilist = &mut self.parse;
        for i in 0..n {
            let bwt_char = self.bwt[i] as usize;
            ilist[frequency_vector[bwt_char] as usize] = i as u32;
            frequency_vector[bwt_char] += 1;
        }

        let mut output = File::create(self.basename.with_extension(".ilist"))
            .map(BufWriter::new)
            .unwrap();
        let result = ilist
            .iter()
            .map(|&i| output.write_u32::<NativeEndian>(i))
            .collect::<Result<(), _>>();
        result.expect("Unable to write ilist vector to file");
    }
}

fn compute_suffix_array_and_lcp(dictionary: &[u8]) -> (Vec<u32>, Vec<i32>) {
    let mut sa = Vec::new();
    let mut lcp = Vec::new();

    sa.resize(dictionary.len(), 0);
    lcp.resize(dictionary.len(), 0);

    unsafe {
        gsacak_u8(dictionary, &mut sa, Some(&mut lcp));
    }

    (sa, lcp)
}

fn get_suffix_length(suffix: u32, suffix_array: &[u32]) -> (u32, usize) {
    let pos = suffix_array
        .binary_search(&suffix)
        .err()
        .expect("suffix was found in suffix array");
    let ret = if pos == suffix_array.len() { (0, 0) } else { (suffix_array[pos] - suffix, pos) };

    ret
}

pub struct PFBwt {
    dictionary: Vec<u8>,
    istart: Vec<u32>,
    ilist: Vec<u32>,
    bwlast: Vec<u8>,
    bwt_file: BufWriter<File>,
}

impl PFBwt {
    pub fn new(basename: &Path) -> PFBwt {
        let mut dictionary = Vec::new();
        read_file_into_vec_u8(basename.with_extension(".dict"), &mut dictionary);
        dictionary.truncate(dictionary.len() - 1);

        let mut istart = Vec::new();
        read_file_into_vec_u32(basename.with_extension(".occ"), &mut istart);

        let mut ilist = Vec::new();
        read_file_into_vec_u32(basename.with_extension(".ilist"), &mut ilist);
        ilist.truncate(ilist.len() - 1);

        let mut bwlast = Vec::new();
        read_file_into_vec_u8(basename.with_extension(".bwlast"), &mut bwlast);
        bwlast.truncate(bwlast.len() - 1);

        let mut last = 1;
        for occ in istart.iter_mut() {
            (*occ, last) = (last, last + *occ);
        }

        let bwt_file = File::create(basename.with_extension(".bwt"))
            .map(BufWriter::new)
            .unwrap();

        PFBwt {
            dictionary,
            istart,
            ilist,
            bwlast,
            bwt_file,
        }
    }

    pub fn contruct_bwt(&mut self, window: u32) {
        let dsize = self.dictionary.len();
        let dwords = self.istart.len() - 1;
        let (dict_sa, dict_lcp) = compute_suffix_array_and_lcp(&self.dictionary);

        self.dictionary[0] = 0;
        let eos = &dict_sa[1..=dwords];

        let mut full_words = 0;
        let mut easy_bwts = 0;
        let mut hard_bwts = 0;

        let mut i = dwords + window as usize + 1;
        while i < dsize {
            let mut next = i + 1;
            let (suffix_len, mut seq_id) = get_suffix_length(dict_sa[i], eos);
            if suffix_len <= window {
                i = next;
                continue;
            }

            if dict_sa[i] == 0 || self.dictionary[dict_sa[i] as usize - 1] == END_OF_WORD {
                full_words += 1;

                for j in self.istart[seq_id]..self.istart[seq_id + 1] {
                    let next_bwt = self.bwlast[self.ilist[j as usize] as usize];
                    self.bwt_file.write(&[next_bwt]).unwrap();
                    easy_bwts += 1;
                }
                i = next;
                continue;
            }

            let mut ids_to_merge = vec![seq_id];
            let mut chars_to_write = vec![self.dictionary[dict_sa[i] as usize - 1]];

            while next < dsize && dict_lcp[next] >= suffix_len as i32 {
                let (next_suffix_len, next_seq_id) = get_suffix_length(dict_sa[next], eos);

                seq_id = next_seq_id;
                if next_suffix_len != suffix_len {
                    break;
                }

                ids_to_merge.push(seq_id);
                chars_to_write.push(self.dictionary[dict_sa[next] as usize - 1]);
                next += 1;
            }

            self.write_chars_with_same_suffix(
                &ids_to_merge,
                &chars_to_write,
                &mut easy_bwts,
                &mut hard_bwts,
            );

            i = next;
        }

        eprintln!("Full words: {full_words}");
        eprintln!("Easy BWT chars: {easy_bwts}");
        eprintln!("Hard BWT chars: {hard_bwts}");
    }

    fn write_chars_with_same_suffix(
        &mut self,
        ids_to_merge: &[usize],
        chars_to_write: &[u8],
        easy_bwts: &mut usize,
        hard_bwts: &mut usize,
    ) {
        let num_words = ids_to_merge.len();
        let mut same_char = true;

        let mut i = 1;
        while i < num_words && same_char {
            same_char = chars_to_write[i - 1] == chars_to_write[i];
            i += 1;
        }
        if same_char {
            let c = chars_to_write[0];
            for &id in ids_to_merge.iter() {
                for _ in self.istart[id]..self.istart[id + 1] {
                    self.bwt_file.write(&[c]).unwrap();
                }
                *easy_bwts += (self.istart[id + 1] - self.istart[id]) as usize;
            }
        } else {
            let mut heap = BinaryHeap::with_capacity(ids_to_merge.len());

            for (i, &id) in ids_to_merge.iter().enumerate() {
                let remaining = self.istart[id + 1] - self.istart[id];
                let bwt_positions = &self.ilist[self.istart[id] as usize..];
                let char_to_write = chars_to_write[i];

                heap.push(Reverse(SeqID {
                    // id,
                    remaining,
                    bwt_positions,
                    char_to_write,
                }));
            }

            while !heap.is_empty() {
                let mut s = heap.pop().unwrap();
                self.bwt_file.write(&[s.0.char_to_write]).unwrap();

                *hard_bwts += 1;
                if s.0.next() {
                    heap.push(s);
                }
            }
        }
    }
}

struct SeqID<'a> {
    // id: usize,
    remaining: u32,
    bwt_positions: &'a [u32],
    char_to_write: u8,
}

impl<'a> SeqID<'a> {
    fn next(&mut self) -> bool {
        self.remaining -= 1;
        self.bwt_positions = &self.bwt_positions[1..];

        self.remaining > 0
    }
}

impl<'a> Eq for SeqID<'a> {}

impl<'a> PartialEq for SeqID<'a> {
    fn eq<'b>(&self, other: &SeqID<'b>) -> bool {
        self.bwt_positions[0] == other.bwt_positions[0]
    }
}

impl<'a> PartialOrd for SeqID<'a> {
    fn partial_cmp<'b>(&self, other: &SeqID<'b>) -> Option<Ordering> {
        let (left, right) = (self.bwt_positions[0], other.bwt_positions[0]);

        left.partial_cmp(&right)
    }
}

impl<'a> Ord for SeqID<'a> {
    fn cmp<'b>(&self, other: &SeqID<'b>) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}
