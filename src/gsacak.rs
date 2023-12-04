use std::collections::HashMap;

const EMPTY: i32 = 1 << 31;

pub type Table = HashMap<u8, usize>;

pub trait Alphabet<T> {
    fn get(&self, index: usize) -> usize;
}

pub fn count_characters(text: &[u8], table: &Table) -> Vec<u32> {
    text.iter().fold(vec![0; table.len()], |mut counts, c| {
        let a = table.get(c).expect("Unknown character in text");
        counts[*a] += 1;

        counts
    })
}

pub fn get_buckets(counts: &[u32], buckets: &mut [u32], end: bool) {
    let mut sum = 0;

    buckets.fill(0);
    for (count, bucket) in counts.iter().zip(buckets.iter_mut()) {
        if end {
            sum += count;
            *bucket = sum - 1;
        } else {
            *bucket = sum;
            sum += count;
        }
    }
}

fn put_lms_substring_level0(
    text: &[u8],
    suffix_array: &mut [u32],
    buckets: &mut [u32],
    table: &Table,
) {
    let n = text.len();
    let mut previous_is_s_type = false;
    for i in (1..(n - 1)).rev() {
        let c1 = table.get(&text[i - 1]).unwrap();
        let c2 = table.get(&text[i]).unwrap();

        let current_is_s_type = *c1 < *c2 || *c1 == *c2 && previous_is_s_type;
        if !current_is_s_type && previous_is_s_type {
            suffix_array[buckets[*c2] as usize] = i as u32;
            buckets[*c2] -= 1;
        }
        previous_is_s_type = current_is_s_type;
    }

    suffix_array[0] = (n as u32) - 1;
}

fn put_lms_substring_leveln(text: &[u32], suffix_array: &mut [u32]) {
    let n = text.len();

    suffix_array.fill(EMPTY as u32);
    let mut current_is_s_type = false;
    let mut current_character = text[n - 2] as usize;

    for i in (1..(n - 1)).rev() {
        let next_character = current_character;
        let successor_is_s_type = current_is_s_type;

        current_character = text[i - 1] as usize;
        current_is_s_type = current_character < next_character
            || (current_character == next_character && successor_is_s_type);

        if !current_is_s_type && successor_is_s_type {
            if suffix_array[next_character] as i32 >= 0 {
                let mut prev = suffix_array[next_character];
                let mut h = next_character + 1;
                let iter = suffix_array[h..]
                    .iter_mut()
                    .take_while(|&&mut s| s as i32 >= 0);
                for s in iter {
                    (*s, prev) = (prev, *s);
                    h += 1;
                }
                suffix_array[h] = prev;
                suffix_array[next_character] = EMPTY as u32;
            }

            let d = suffix_array[next_character] as i32;
            if d == EMPTY {
                if suffix_array[next_character - 1] as i32 == EMPTY {
                    suffix_array[next_character] = -1i32 as u32;
                    suffix_array[next_character - 1] = i as u32;
                } else {
                    suffix_array[next_character] = i as u32;
                }
            } else {
                let mut position = next_character - d.unsigned_abs() as usize - 1;
                if suffix_array[position] as i32 != EMPTY {
                    position += 1;
                    suffix_array[position..=next_character].rotate_right(1);
                } else {
                    suffix_array[next_character] -= 1;
                }
                suffix_array[position] = i as u32;
            }
        }
    }

    for i in (1..n).rev() {
        let j = suffix_array[i] as i32;
        if j < 0 && j != EMPTY {
            let start = i - j.unsigned_abs() as usize;
            suffix_array[start..=i].rotate_right(1);
            suffix_array[start] = EMPTY as u32;
        }
    }

    suffix_array[0] = n as u32 - 1;
}

fn put_suffix_level0(text: &[u8], suffix_array: &mut [u32], buckets: &mut [u32], table: &Table) {
    let n = suffix_array.len();

    for i in (1..n).rev() {
        let j = suffix_array[i] as usize;
        suffix_array[i] = 0;
        let c = table.get(&text[j]).unwrap();
        suffix_array[buckets[*c] as usize] = j as u32;
        buckets[*c] -= 1;
    }
    suffix_array[0] = text.len() as u32 - 1;
}

fn induce_l_type_level0(
    text: &[u8],
    suffix_array: &mut [u32],
    buckets: &mut [u32],
    table: &Table,
    suffix: bool,
) {
    let n = text.len() - 1;

    for i in 0..n {
        if suffix_array[i] == 0 {
            continue;
        }
        let j = (suffix_array[i] - 1) as usize;
        let c1 = table.get(&text[j]).unwrap();
        let c2 = table.get(&text[j + 1]).unwrap();
        if *c1 >= *c2 {
            suffix_array[buckets[*c1] as usize] = j as u32;
            buckets[*c1] += 1;
            if !suffix && i > 0 {
                suffix_array[i] = 0;
            }
        }
    }
}

fn induce_s_type_level0(
    text: &[u8],
    suffix_array: &mut [u32],
    buckets: &mut [u32],
    table: &Table,
    suffix: bool,
) {
    let n = text.len();

    for i in (1..n).rev() {
        if suffix_array[i] == 0 {
            continue;
        }
        let j = (suffix_array[i] - 1) as usize;
        let c1 = table.get(&text[j]).unwrap();
        let c2 = table.get(&text[j + 1]).unwrap();
        if *c1 <= *c2 && (buckets[*c1] as usize) < i {
            suffix_array[buckets[*c1] as usize] = j as u32;
            buckets[*c1] -= 1;
            if !suffix {
                suffix_array[i] = 0;
            }
        }
    }
}

fn induce_l_type_leveln(text: &[u32], suffix_array: &mut [u32], suffix: bool) {
    let n = text.len();
    let mut step;
    let mut i = 0;

    while i < n {
        step = 1;
        if (suffix_array[i] as i32) <= 0 {
            i += step;
            continue;
        }
        let j = (suffix_array[i] as i32 - 1) as usize;
        let (current, next) = (text[j] as usize, text[j + 1] as usize);
        if current < next {
            i += step;
            continue;
        }

        let mut d = suffix_array[current] as i32;
        if d >= 0 {
            let mut prev = suffix_array[current];
            let mut counter = current - 1;
            let bucket_iter = suffix_array[0..current]
                .iter_mut()
                .rev()
                .take_while(|&&mut s| (s as i32) >= 0 || s as i32 == EMPTY);
            for s in bucket_iter {
                (prev, *s) = (*s, prev);
                counter -= 1;
            }

            suffix_array[counter] = prev;
            d = EMPTY;
            if counter < i {
                step = 0;
            }
        }

        if d == EMPTY {
            if current < n - 1 && suffix_array[current + 1] == EMPTY as u32 {
                suffix_array[current] = -1i32 as u32;
                suffix_array[current + 1] = j as u32;
            } else {
                suffix_array[current] = j as u32;
            }
        } else {
            let mut position = current + d.unsigned_abs() as usize + 1;
            if position > n - 1 || (suffix_array[position] as i32) != EMPTY {
                position -= 1;
                suffix_array[current..=position].rotate_left(1);
                if current < i {
                    step = 0;
                }
            } else {
                suffix_array[current] -= 1;
            }
            suffix_array[position] = j as u32;
        }
        let l_type = (j + 1 < n - 1)
            && ((text[j + 1] > text[j + 2])
                || (text[j + 1] == text[j + 2] && (text[j + 1] as usize) < i));
        if (!suffix || !l_type) && i > 0 {
            let k = if step == 0 { i - 1 } else { i };
            suffix_array[k] = EMPTY as u32;
        }
        i += step;
    }
    for i in 1..n {
        let j = suffix_array[i] as i32;
        if j < 0 && j != EMPTY {
            let start = i;
            let stop = i + j.unsigned_abs() as usize;
            suffix_array[start..=stop].rotate_left(1);
            suffix_array[stop] = EMPTY as u32;
        }
    }
}

fn induce_s_type_leveln(text: &[u32], suffix_array: &mut [u32], suffix: bool) {
    let n = text.len();
    let mut i = n - 1;
    let mut step;

    while i > 0 {
        step = 1;
        if (suffix_array[i] as i32) <= 0 {
            i -= step;
            continue;
        }
        let j = suffix_array[i] as usize - 1;
        let (current, next) = (text[j] as usize, text[j + 1] as usize);
        let is_s_type = current < next || (current == next && current > i);
        if !is_s_type {
            i -= step;
            continue;
        }

        let mut d = suffix_array[current] as i32;

        if d >= 0 {
            let mut prev = suffix_array[current];
            let mut counter = current + 1;
            let iter = &mut suffix_array[counter..]
                .iter_mut()
                .take_while(|&&mut s| (s as i32) >= 0 || s as i32 == EMPTY);
            for s in iter {
                (*s, prev) = (prev, *s);
                counter += 1;
            }
            suffix_array[counter] = prev;
            if counter > i {
                step = 0;
            }
            d = EMPTY;
        }

        if d == EMPTY {
            if suffix_array[current - 1] as i32 == EMPTY {
                suffix_array[current] = -1i32 as u32;
                suffix_array[current - 1] = j as u32;
            } else {
                suffix_array[current] = j as u32;
            }
        } else {
            let mut position = current - d.unsigned_abs() as usize - 1;
            if suffix_array[position] as i32 != EMPTY {
                position += 1;
                suffix_array[position..=current].rotate_right(1);
                if current > i {
                    step = 0;
                }
            } else {
                suffix_array[current] -= 1;
            }
            suffix_array[position] = j as u32;
        }
        if !suffix {
            let j = if step == 0 { i + 1 } else { i };
            suffix_array[j] = EMPTY as u32;
        }
        i -= step;
    }

    for i in (1..(suffix_array.len() - 1)).rev() {
        let j = suffix_array[i] as i32;
        if j < 0 && j != EMPTY {
            let start = i - j.unsigned_abs() as usize;
            let stop = i;
            suffix_array[start..=stop].rotate_right(1);
            suffix_array[start] = EMPTY as u32;
        }
    }
}

fn name_lms_substrings<T: PartialOrd + Debug>(text: &[T], suffix_array: &mut [u32], n1: usize) -> usize {
    let mut name_counter = 0;

    {
        let (suffix_array1, name_buffer) = suffix_array.split_at_mut(n1);
        name_buffer.fill(EMPTY as u32);

        let mut name: usize = 0;
        let mut pre_position = 0;
        let mut pre_length = 0;

        for i in 0..n1 {
            let position = suffix_array1[i] as usize;
            let length = get_length_of_lms_substring(&text[position..]);
            let diff = length != pre_length
                || text[position..position + length] != text[pre_position..pre_position + length];

            if diff {
                name = i;
                name_counter += 1;
                suffix_array1[name] = 1;
                (pre_position, pre_length) = (position, length);
            } else {
                suffix_array1[name] += 1;
            }
            name_buffer[position / 2] = name as u32;
        }

        let mut j = name_buffer.len() - 1;
        for i in (0..name_buffer.len()).rev() {
            if name_buffer[i] as i32 != EMPTY {
                name_buffer[j] = name_buffer[i];
                j -= 1;
            }
        }
    }

    let mut successor_is_s_type = true;
    let (sa1, t1) = suffix_array.split_at_mut(suffix_array.len() - n1);

    for i in (1..n1).rev() {
        let successor = t1[i];
        let current = t1[i - 1];
        let current_is_s_type =
            current < successor || (current == successor && successor_is_s_type);
        if current_is_s_type {
            t1[i - 1] += sa1[t1[i - 1] as usize] - 1;
        }
        successor_is_s_type = current_is_s_type;
    }

    name_counter
}

use std::fmt::Debug;
fn get_suffix_array_lms<T: PartialOrd + Debug>(
    text: &[T],
    suffix_array: &mut [u32],
    text1: &mut [u32],
    level0: bool,
) {
    let n1 = text1.len();
    let mut j = n1 - 1;
    let n = text.len();

    text1[j] = (n - 1) as u32;
    j -= 1;
    let mut successor_is_s_type = false;
    for i in (1..(n - 1)).rev() {
        let current_is_s_type =
            text[i - 1] < text[i] || (text[i - 1] == text[i] && successor_is_s_type);
        if !current_is_s_type && successor_is_s_type {
            text1[j] = i as u32;
            j = j.saturating_sub(1);
        }
        successor_is_s_type = current_is_s_type;
    }
    for i in 0..n1 {
        suffix_array[i] = text1[suffix_array[i] as usize];
    }

    suffix_array[n1..].fill(if level0 { 0 } else { EMPTY as u32 });
    text1.fill(if level0 { 0 } else { EMPTY as u32 });
}

fn get_length_of_lms_substring<T: PartialOrd>(text: &[T]) -> usize {
    let n = text.len();

    if n == 1 {
        return 1;
    }

    let mut lms_length = 0;
    let mut i = 1;

    loop {
        if text[i] < text[i - 1] {
            break;
        }
        i += 1;
    }

    loop {
        if (i > n - 1) || (text[i] > text[i - 1]) {
            break;
        }
        if (i == n - 1) || (text[i] < text[i - 1]) {
            lms_length = i;
        }
        i += 1;
    }

    lms_length + 1
}

fn sacak_leveln(text: &[u32], suffix_array: &mut [u32]) {
    put_lms_substring_leveln(text, suffix_array);
    induce_l_type_leveln(text, suffix_array, false);
    induce_s_type_leveln(text, suffix_array, false);
    dbg!(&suffix_array);
    let mut n1 = 0;
    for i in 0..suffix_array.len() {
        if suffix_array[i] as i32 > 0 {
            suffix_array[n1] = suffix_array[i];
            n1 += 1;
        }
    }

    let name_counter = name_lms_substrings(text, suffix_array, n1);
    let (suffix_array1, text1) = suffix_array.split_at_mut(suffix_array.len() - n1);

    if name_counter < n1 {
        sacak_leveln(text1, suffix_array1);
    } else {
        for i in 0..n1 {
            suffix_array1[text1[i] as usize] = i as u32;
        }
    }
    get_suffix_array_lms(text, suffix_array1, text1, false);
    put_suffix_leveln(text, suffix_array, n1);
    induce_l_type_leveln(text, suffix_array, true);
    induce_s_type_leveln(text, suffix_array, true);
}

fn put_suffix_leveln(text: &[u32], suffix_array: &mut [u32], n1: usize) {
    let mut prev = -1;
    let mut pos = 0;

    // dbg!(&suffix_array);
    // dbg!(&text);
    for i in (1..n1).rev() {
        let j = suffix_array[i];

        suffix_array[i] = EMPTY as u32;
        let curr = text[j as usize] as i32;

        if curr != prev {
            prev = curr;
            pos = curr as usize;
        }
        suffix_array[pos] = j;
        pos -= 1;
    }
}

fn sacak0(text: &[u8], suffix_array: &mut [u32], table: Table) {
    let text = text.as_ref();
    let counts = count_characters(text, &table);
    let mut buckets = vec![0; counts.len()];

    get_buckets(&counts, &mut buckets, true);
    put_lms_substring_level0(text, suffix_array, &mut buckets, &table);

    get_buckets(&counts, &mut buckets, false);
    induce_l_type_level0(text, suffix_array, &mut buckets, &table, false);

    get_buckets(&counts, &mut buckets, true);
    induce_s_type_level0(text, suffix_array, &mut buckets, &table, false);

    let mut n1 = 0;
    for i in 0..text.len() {
        if suffix_array[i] > 0 {
            suffix_array[n1] = suffix_array[i];
            n1 += 1;
        }
    }

    let name_counter = name_lms_substrings(text, suffix_array, n1);
    let (suffix_array1, text1) = suffix_array.split_at_mut(suffix_array.len() - n1);

    if name_counter < n1 {
        sacak_leveln(text1, suffix_array1);
    } else {
        for i in 0..n1 {
            suffix_array1[text1[i] as usize] = i as u32;
        }
    }
    get_suffix_array_lms(text, suffix_array1, text1, true);

    get_buckets(&counts, &mut buckets, true);
    put_suffix_level0(text, &mut suffix_array[0..n1], &mut buckets, &table);

    get_buckets(&counts, &mut buckets, false);
    induce_l_type_level0(text, suffix_array, &mut buckets, &table, true);

    get_buckets(&counts, &mut buckets, true);
    induce_s_type_level0(text, suffix_array, &mut buckets, &table, true);
}

pub fn sacak<T: AsRef<[u8]>>(text: T, table: Table) -> Vec<u32> {
    let text = text.as_ref();

    if text.len() <= 1 {
        Vec::new()
    } else {
        let mut suffix_array = vec![0; text.len()];
        sacak0(text, &mut suffix_array, table);

        suffix_array
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_chars() {
        // let map: Vec<(u8, usize)> = vec![
        //     ('\0' as u8, 0),
        //     ('A' as u8, 1),
        //     ('C' as u8, 2),
        //     ('G' as u8, 3),
        //     ('T' as u8, 4),
        // ];
        let map: Vec<(u8, usize)> = vec![
            ('\0' as u8, 0),
            ('I' as u8, 1),
            ('M' as u8, 2),
            ('P' as u8, 3),
            ('S' as u8, 4),
        ];

        let table = map.into_iter().collect::<Table>();
        // let data = "TTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGGGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCATCTGGTAGCGATGATTGA\0";
        let data = "MMIISSIISSIIPPII\0";

        let counts = count_characters(data.as_bytes(), &table);
        let sa = sacak(data, table);
        assert_eq!(counts[0], 17);
    }

    // #[test]
    // fn test_lms_length() {
    //     let string = "cinnut\0";

    //     assert_eq!(3, get_length_of_lms_substring(&string.as_bytes()));
    // }
}
