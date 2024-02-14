use std::fmt::Debug;

trait Character: Copy + Ord + Debug {
    fn to_usize(self) -> usize;

    fn zero() -> Self;
}

impl Character for i32 {
    #[inline]
    fn to_usize(self) -> usize {
        self as usize
    }

    #[inline]
    fn zero() -> i32 {
        0
    }
}

impl Character for u8 {
    #[inline]
    fn to_usize(self) -> usize {
        self as usize
    }

    #[inline]
    fn zero() -> u8 {
        0
    }
}

#[inline]
fn get_buckets<C: Character>(text: &[C], buckets: &mut [i32], characters: usize, end: bool) {
    buckets.fill(0);
    assert!(buckets.len() == characters);

    for &c in text {
        buckets[c.to_usize()] += 1;
    }

    let mut sum = 0;

    if end {
        for c in buckets {
            sum += *c;
            *c = sum;
        }
    } else {
        for c in buckets {
            sum += *c;
            *c = sum - *c;
        }
    }
}

#[inline]
fn induce_sa<C: Character + Debug>(
    text: &[C],
    suffix_array: &mut [i32],
    buckets: &mut [i32],
    alpha_size: usize,
    suffix: bool
) {
    get_buckets(text, buckets, alpha_size, false);

    let n = text.len();
    let mut j = n as i32 - 1;
    let mut c1 = text[j as usize];
    let mut pos = buckets[c1.to_usize()] as usize;

    suffix_array[pos] = if j > 0 && text[j as usize - 1] < c1 {
        !j
    } else {
        j
    };
    pos += 1;

    for i in 0..n {
        j = suffix_array[i];
        suffix_array[i] = !j;

        if j > 0 {
            j -= 1;
            let c0 = text[j as usize];

            if c0 != c1 {
                buckets[c1.to_usize()] = pos as i32;
                c1 = c0;
                pos = buckets[c1.to_usize()] as usize;
            }

            suffix_array[pos] = if j > 0 && text[j as usize - 1] < c1 {
                !j
            } else {
                j
            };
            pos += 1;
        }
    }

    get_buckets(text, buckets, alpha_size, true);
    c1 = C::zero();
    pos = buckets[c1.to_usize()] as usize;

    for i in (0..n).rev() {
        let mut j = suffix_array[i];

        if j > 0 {
            j -= 1;
            let c0 = text[j as usize];
            if c0 != c1 {
                buckets[c1.to_usize()] = pos as i32;
                c1 = c0;
                pos = buckets[c1.to_usize()] as usize;
            }
            pos -= 1;
            suffix_array[pos] = if j == 0 || text[j as usize - 1] > c1 {
                !j
            } else {
                j
            };
        } else {
            suffix_array[i] = !j;
        }
    }
}

fn sais<C: Character>(
    text: &[C],
    mut suffix_array: &mut [i32],
    alphabet_size: usize,
    free_space: usize,
) {
    let buckets;
    if alphabet_size <= free_space {
        (suffix_array, buckets) = suffix_array.split_at_mut(suffix_array.len() - alphabet_size);
    } else {
        buckets = unsafe {
            let ptr = vec![0i32; alphabet_size].as_mut_ptr();
            std::slice::from_raw_parts_mut(ptr, alphabet_size)
        }
    }

    get_buckets(text, buckets, alphabet_size, true);
    let n = text.len();
    suffix_array.fill(0);

    // place the LMS suffixes
    let mut c1 = text[n - 1];
    let mut s_type = false;

    for i in (0..n - 1).rev() {
        let c0 = text[i];
        if c0 < c1 {
            s_type = true;
        } else if c0 > c1 && s_type {
            buckets[c1.to_usize()] -= 1;
            suffix_array[buckets[c1.to_usize()] as usize] = (i + 1) as i32;
            s_type = false;
        }
        c1 = c0;
    }

    // induce S and L type suffixes
    induce_sa(text, suffix_array, buckets, alphabet_size, false);

    // Store the sorted LMS-substrings in the first m items of SA
    let mut m: usize = 0;
    for i in 0..n {
        let pos = suffix_array[i];
        let c0 = text[pos as usize];

        if pos > 0 && text[pos as usize - 1] > c0 {
            let mut j = pos as usize + 1;

            while j < n {
                c1 = text[j];

                if c0 != c1 {
                    break;
                }

                j += 1;
            }

            if j < n && c0 < c1 {
                suffix_array[m] = pos;
                m += 1;
            }
        }
    }

    // for i in 0..n {
    //     let j = suffix_array[i];
    //     suffix_array[i] = 0;
    //     if j != 0 {
    //         suffix_array[m] = j;
    //         m += 1;
    //     }
    // }

    let mut name = 0;
    {
        let (sa, name_buffer) = suffix_array.split_at_mut(m);
        name_buffer.fill(0);

        let mut j = n;
        let mut c1 = text[j - 1];
        let mut is_stype = false;

        for i in (0..n - 1).rev() {
            let c0 = text[i];
            if c0 < c1 {
                is_stype = true;
            } else if c0 > c1 && is_stype {
                name_buffer[(i + 1) >> 1] = (j - i - 1) as i32;
                is_stype = false;
                j = i + 1;
            }
            c1 = c0;
        }

        let mut prev_pos = n;
        let mut prev_len = 0;

        for i in 0..m {
            let curr_pos = sa[i] as usize;
            let curr_len = name_buffer[curr_pos >> 1] as usize;
            let mut diff = true;

            if prev_len == curr_len {
                diff = text[curr_pos..curr_pos + curr_len] != text[prev_pos..prev_pos + prev_len];
            }

            if diff {
                name += 1;
                prev_pos = curr_pos;
                prev_len = curr_len;
            }

            name_buffer[curr_pos >> 1] = name;
        }
    }

    if (name as usize) < m {
        let mut j = suffix_array.len() - 1;

        for i in (m..n).rev() {
            if suffix_array[i] != 0 {
                suffix_array[j] = suffix_array[i] - 1;
                j = j.saturating_sub(1);
            }
        }

        let (suffix_array, reduced_text) = suffix_array.split_at_mut(suffix_array.len() - m);
        let free_space = free_space + n - m * 2;

        sais(reduced_text, suffix_array, name as usize, free_space);

        let mut j = m;
        let mut c1 = text[n - 1];
        let mut s_type = false;

        for i in (0..n - 1).rev() {
            let c0 = text[i];
            if c0 < c1 {
                s_type = true;
            } else if c0 > c1 && s_type {
                reduced_text[j-1] = i as i32 + 1;
                j -= 1;
                s_type = false;
            }
            c1 = c0;
        }

        for i in 0..m {
            suffix_array[i] = reduced_text[suffix_array[i] as usize];
        }
    }

    get_buckets(text, buckets, alphabet_size, true);
    suffix_array[m..].fill(0);
    for i in (0..m).rev() {
        let j = suffix_array[i];
        let c = text[j as usize].to_usize();
        suffix_array[i] = 0;
        buckets[c] -= 1;
        suffix_array[buckets[c] as usize] = j;
    }

    induce_sa(text, suffix_array, buckets, alphabet_size, true);
}

pub fn construct_suffix_array(text: impl AsRef<[u8]>, alphabet_size: usize) -> Vec<i32> {
    let text = text.as_ref();
    let n = text.len();

    if n == 0 {
        return vec![];
    }

    let mut suffix_array = vec![0; n + 1];
    suffix_array[0] = n as i32;

    if n > 1 {
        sais(text, &mut suffix_array[1..], alphabet_size, 0);
    }

    suffix_array
}

#[cfg(test)]
mod tests {
    use super::construct_suffix_array;

    #[test]
    fn try_suffix_array() {
        let text = "TTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGGGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCATCTGGTAGCGATGATTGA";
        let sa = construct_suffix_array(text, 256);
        let result = &[
            140, 139, 108, 58, 109, 20, 117, 28, 59, 110, 21, 63, 9, 114, 118, 87, 29, 26, 60, 128,
            69, 73, 79, 14, 49, 111, 11, 121, 132, 22, 64, 51, 135, 107, 19, 116, 62, 113, 86, 25,
            68, 72, 78, 13, 10, 120, 106, 115, 85, 67, 119, 100, 101, 102, 46, 103, 88, 92, 47,
            130, 104, 2, 39, 5, 89, 30, 93, 95, 97, 36, 54, 123, 138, 57, 27, 48, 131, 134, 18, 61,
            24, 71, 77, 105, 84, 66, 99, 45, 129, 38, 4, 56, 17, 70, 76, 83, 3, 75, 74, 40, 6, 125,
            80, 41, 7, 126, 90, 15, 81, 42, 31, 8, 127, 50, 112, 12, 91, 1, 94, 96, 35, 53, 122,
            137, 133, 23, 65, 98, 44, 37, 55, 16, 82, 124, 0, 34, 52, 136, 43, 33, 32,
        ];

        assert_eq!(&sa, result);
    }
}
