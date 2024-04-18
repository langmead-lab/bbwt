use crate::utils::compute_lcp_phi_sparse;
use crate::int::Character;
use crate::utils::StackRMQ;
use crate::int::Index;

#[inline]
fn get_buckets<C: Character, I: Index>(text: &[C], buckets: &mut [I], characters: usize, end: bool) {
    buckets.fill(I::ZERO);
    assert!(buckets.len() == characters);

    for &c in text {
        buckets[c.char_to_usize()] += I::ONE;
    }

    let mut sum = I::ZERO;

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

fn induce_sa_and_lcp_generalized<C: Character, I: Index>(
    text: &[C],
    suffix_array: &mut [I],
    lcp: &mut [I],
    buckets: &mut [I],
    separator: C,
    alpha_size: usize,
) {
    let n = suffix_array.len();
    // mark the positions of the L/S seam in each bucket
    for i in separator.char_to_usize() + 1..alpha_size {
        let pos = buckets[i].as_usize();
        if pos < n && suffix_array[pos] != I::max_value() {
            lcp[pos] = I::max_value();
        }
    }

    get_buckets(text, buckets, alpha_size, false);

    // Document separator is at the end of the string, we no
    // longer need to do this since the get skipped anyway
    let mut j = I::from(n).unwrap() - I::ONE;
    let mut c1 = text[j.as_usize()];
    let mut pos = buckets[c1.char_to_usize()].as_usize();
    let mut rmq = StackRMQ::new(alpha_size, false);

    // suffix_array[pos] = if j > 0 && text[j as usize - 1] < c1 {
    //     !j
    // } else {
    //     j
    // };
    // pos += 1;

    for i in 0..n {
        if suffix_array[i] == I::max_value() {
            continue;
        }

        j = suffix_array[i];
        suffix_array[i] = !j;

        if lcp[i] == I::max_value() {
            // we've encountered an L/S seam calculate LCP by direct comparison
            let s1 = if suffix_array[i] < I::ZERO {
                !suffix_array[i]
            } else {
                suffix_array[i]
            };
            // find the position of the last L-type suffix placed in this c-bucket.
            let c0 = text[s1.as_usize()];
            let k = if c0 == c1 {
                pos - 1
            } else {
                buckets[c0.char_to_usize()].as_usize() - 1
            };
            let s2 = if suffix_array[k] < I::ZERO {
                !suffix_array[k]
            } else {
                suffix_array[k]
            };
            let iter1 = text[s1.as_usize()..].iter();
            let iter2 = text[s2.as_usize()..].iter();

            lcp[i] = I::from(iter1.zip(iter2).take_while(|(&a, &b)| a == b).count()).unwrap();
            lcp[k] = !lcp[k];
        }

        rmq.push(i, lcp[i]);

        if j > I::ZERO {
            j -= I::ONE;

            let c0 = text[j.as_usize()];

            if c0 == separator {
                continue;
            }

            if c0 != c1 {
                buckets[c1.char_to_usize()] = I::from(pos).unwrap();
                c1 = c0;
                pos = buckets[c1.char_to_usize()].as_usize();
            }

            suffix_array[pos] = if j > I::ZERO && text[j.as_usize() - 1] < c1 {
                !j
            } else {
                j
            };

            let min_lcp = rmq.get_min(c0, i);
            lcp[pos] = min_lcp;
            pos += 1;
        }
    }

    get_buckets(text, buckets, alpha_size, true);
    c1 = C::zero();
    pos = buckets[c1.char_to_usize()].as_usize();
    rmq.reset_for_type(true);

    for i in (0..n).rev() {
        let mut j = suffix_array[i];

        if j > I::ZERO {
            j -= I::ONE;

            let c0 = text[j.as_usize()];

            if c0 == separator {
                continue;
            }

            let min_lcp = rmq.get_min(c0, i);

            if c0 != c1 {
                buckets[c1.char_to_usize()] = I::from(pos).unwrap();
                c1 = c0;
                pos = buckets[c1.char_to_usize()].as_usize();
            }
            pos -= 1;
            suffix_array[pos] = if j == I::ZERO || text[j.as_usize() - 1] > c1 {
                !j
            } else {
                j
            };

            lcp[pos + 1] = min_lcp;
        } else {
            suffix_array[i] = !j;
        }

        if i > 0 && lcp[i - 1] < I::ZERO {
            // We have encountered an L/S seam when placing the first S-type suffix.
            // Calculate the LCP by direct comparison.
            let s1 = if suffix_array[i] < I::ZERO {
                !suffix_array[i].as_usize()
            } else {
                suffix_array[i].as_usize()
            };
            let s2 = if suffix_array[i - 1] < I::ZERO {
                (!suffix_array[i - 1]).as_usize()
            } else {
                suffix_array[i - 1].as_usize()
            };
            let iter1 = text[s1..].iter();
            let iter2 = text[s2..].iter();

            lcp[i] = I::from(iter1.zip(iter2).take_while(|(&a, &b)| a == b).count()).unwrap();
            lcp[i - 1] = !lcp[i - 1];
        }

        rmq.push(i, lcp[i]);
    }
}

fn induce_sa_generalized<C: Character, I: Index>(
    text: &[C],
    suffix_array: &mut [I],
    buckets: &mut [I],
    separator: C,
    alpha_size: usize,
    induce_lms: bool,
) {
    get_buckets(text, buckets, alpha_size, false);

    let n = suffix_array.len();
    // Document separator is at the end of the string, we no
    // longer need to do this since the get skipped anyway
    let mut j = I::from(n - 1).unwrap();
    let mut c1 = text[j.as_usize()];
    let mut pos = buckets[separator.char_to_usize()].as_usize();

    // suffix_array[pos] = if j > 0 && text[j as usize - 1] < c1 {
    //     !j
    // } else {
    //     j
    // };
    // pos += 1;
    // let mut k = 0;
    // for i in 0..n {
    //     let j = suffix_array[i];
    //     // suffix_array[i] = !j;
    //     c1 = text[j as usize];
    //     if c1 != separator {
    //         break;
    //     }
    //     let c0 = text[j as usize - 1];
    //     suffix_array[buckets[c0.char_to_usize()] as usize] = j - 1;
    //     buckets[c0.char_to_usize()] += 1;
    //     k = 3;
    // }

    // get_buckets(text, buckets, alpha_size, false);

    for i in 0..n {
        j = suffix_array[i];
        suffix_array[i] = !j;

        if j > I::ZERO {
            j -= I::ONE;

            let c0 = text[j.as_usize()];

            if c0 == separator {
                continue;
            }

            if c0 != c1 {
                buckets[c1.char_to_usize()] = I::from(pos).unwrap();
                c1 = c0;
                pos = buckets[c1.char_to_usize()].as_usize();
            }

            suffix_array[pos] = if j > I::ZERO && text[j.as_usize() - 1] < c1 {
                !j
            } else {
                j
            };
            pos += 1;
        }
    }

    get_buckets(text, buckets, alpha_size, true);
    c1 = C::zero();
    pos = buckets[c1.char_to_usize()].as_usize();

    for i in (0..n).rev() {
        let mut j = suffix_array[i];

        if j > I::ZERO {
            j -= I::ONE;
            let c0 = text[j.as_usize()];

            if c0 == separator {
                continue;
            }

            if c0 != c1 {
                buckets[c1.char_to_usize()] = I::from(pos).unwrap();
                c1 = c0;
                pos = buckets[c1.char_to_usize()].as_usize();
            }
            pos -= 1;
            suffix_array[pos] = if j == I::ZERO || text[j.as_usize() - 1] > c1 {
                !j
            } else {
                j
            };
        } else {
            suffix_array[i] = !j;
        }
    }
}

#[inline]
fn induce_sa<C: Character, I: Index>(
    text: &[C],
    suffix_array: &mut [I],
    buckets: &mut [I],
    alpha_size: usize,
) {
    get_buckets(text, buckets, alpha_size, false);

    let n = text.len();
    let mut j = I::from(n - 1).unwrap();
    let mut c1 = text[j.as_usize()];
    let mut pos = buckets[c1.char_to_usize()].as_usize();

    suffix_array[pos] = if j > I::ZERO && text[j.as_usize() - 1] < c1 {
        !j
    } else {
        j
    };
    pos += 1;

    for i in 0..n {
        j = suffix_array[i];
        suffix_array[i] = !j;

        if j > I::ZERO {
            j -= I::ONE;
            let c0 = text[j.as_usize()];

            if c0 != c1 {
                buckets[c1.char_to_usize()] = I::from(pos).unwrap();
                c1 = c0;
                pos = buckets[c1.char_to_usize()].as_usize();
            }

            suffix_array[pos] = if j > I::ZERO && text[j.as_usize() - 1] < c1 {
                !j
            } else {
                j
            };
            pos += 1;
        }
    }

    get_buckets(text, buckets, alpha_size, true);
    c1 = C::zero();
    pos = buckets[c1.char_to_usize()].as_usize();

    for i in (0..n).rev() {
        let mut j = suffix_array[i];

        if j > I::ZERO {
            j -= I::ONE;
            let c0 = text[j.as_usize()];
            if c0 != c1 {
                buckets[c1.char_to_usize()] = I::from(pos).unwrap();
                c1 = c0;
                pos = buckets[c1.char_to_usize()].as_usize();
            }
            pos -= 1;
            suffix_array[pos] = if j == I::ZERO || text[j.as_usize() - 1] > c1 {
                !j
            } else {
                j
            };
        } else {
            suffix_array[i] = !j;
        }
    }
}

fn sort_lms_substrings_generalized<C: Character, I: Index>(
    text: &[C],
    suffix_array: &mut [I],
    separator: C,
    buckets: &mut [I],
    alphabet_size: usize,
) -> usize {
    get_buckets(text, buckets, alphabet_size, true);
    let n = text.len();
    suffix_array.fill(I::ZERO);

    buckets[separator.char_to_usize()] -= I::ONE;
    suffix_array[buckets[separator.char_to_usize()].as_usize()] = I::from(n - 1).unwrap();

    // place the LMS suffixes in the end of their buckets
    let mut c1 = text[n - 1];
    let mut s_type = false;
    let mut prev_lms_pos = n - 1;

    for i in (0..n - 1).rev() {
        let c0 = text[i];
        if c0 < c1 {
            s_type = true;
        } else if c0 > c1 && s_type {
            // remove the LMS suffix that induces the separator
            if c1 == separator {
                suffix_array[buckets[text[prev_lms_pos].char_to_usize()].as_usize()] = I::ZERO;
                buckets[text[prev_lms_pos].char_to_usize()] += I::ONE;
            }
            buckets[c1.char_to_usize()] -= I::ONE;
            suffix_array[buckets[c1.char_to_usize()].as_usize()] = I::from(i + 1).unwrap();
            prev_lms_pos = i + 1;
            s_type = false;
        }
        c1 = c0;
    }
    // induce S and L type suffixes
    induce_sa_generalized(text, suffix_array, buckets, separator, alphabet_size, true);

    let mut pos = 0;
    for (i, &c) in text.iter().enumerate() {
        if c == separator {
            suffix_array[pos] = I::from(i).unwrap();
            pos += 1;
        }
    }

    // Store the sorted LMS-substrings in the first m items of SA
    let mut lms_substring_count: usize = 0;
    for i in 0..n {
        let pos = suffix_array[i];
        let c0 = text[pos.as_usize()];

        if pos > I::ZERO {
            // if c0 == separator {
            //     lms_substring_count +=1;
            //     continue;
            // }
            if text[pos.as_usize() - 1] > c0 {
                let mut j = pos.as_usize() + 1;

                while j < n {
                    c1 = text[j];

                    if c0 != c1 {
                        break;
                    }

                    j += 1;
                }

                if j < n && c0 < c1 {
                    suffix_array[lms_substring_count] = pos;
                    lms_substring_count += 1;
                }
            }
        }
    }

    lms_substring_count
}

fn name_lms_substrings_generalized<C: Character, I: Index>(
    text: &[C],
    suffix_array: &mut [I],
    num_lms_substrings: usize,
    lcp: &mut [I],
    separator: C,
) -> usize {
    let n = text.len();
    let mut name_count = 0;
    {
        let (sa, name_buffer) = suffix_array.split_at_mut(num_lms_substrings);
        name_buffer.fill(I::ZERO);

        let mut j = n - 1;
        let mut c1 = text[j];
        let mut is_stype = false;

        for i in (0..n - 1).rev() {
            let c0 = text[i];
            if c0 < c1 {
                is_stype = true;
            } else if c0 > c1 && is_stype {
                name_buffer[(i + 1) >> 1] = I::from(j - i).unwrap();
                is_stype = false;
                j = i + 1;
            }
            c1 = c0;
        }

        let mut prev_pos = n;
        let mut prev_len = 0;

        for i in 0..num_lms_substrings {
            let curr_pos = sa[i].as_usize();
            let curr_len = name_buffer[curr_pos >> 1].as_usize();
            let mut diff = false;
            let iter1 = text[curr_pos..curr_pos + curr_len].iter();
            let iter2 = text[prev_pos..prev_pos + prev_len].iter();
            let curr_lcp = iter1
                .zip(iter2)
                .take_while(|(&a, &b)| a == b && (a != separator || b != separator))
                .count();

            if prev_len != curr_len || curr_lcp < curr_len {
                diff = true;
            }

            if diff {
                name_count += 1;
                prev_pos = curr_pos;
                prev_len = curr_len;
            }

            if lcp.len() > 0 {
                lcp[i] = I::from(curr_lcp).unwrap();
            }
            name_buffer[curr_pos >> 1] = I::from(name_count).unwrap();
        }
    }

    name_count
}

pub fn gsais<C: Character, I: Index>(
    text: &[C],
    mut suffix_array: &mut [I],
    separator: C,
    alphabet_size: usize,
    free_space: usize,
) {
    let n = text.len();
    let buckets;
    if alphabet_size <= free_space {
        (suffix_array, buckets) = suffix_array.split_at_mut(suffix_array.len() - alphabet_size);
    } else {
        buckets = unsafe {
            let ptr = vec![I::ZERO; alphabet_size].as_mut_ptr();
            std::slice::from_raw_parts_mut(ptr, alphabet_size)
        }
    }

    let num_lms_substrings =
        sort_lms_substrings_generalized(text, suffix_array, separator, buckets, alphabet_size);
    let name_count =
        name_lms_substrings_generalized(text, suffix_array, num_lms_substrings, &mut [], separator);

    if (name_count as usize) < num_lms_substrings {
        let mut j = suffix_array.len() - 1;

        for i in (num_lms_substrings..n).rev() {
            if suffix_array[i] != I::ZERO {
                suffix_array[j] = suffix_array[i] - I::ONE;
                j = j.saturating_sub(1);
            }
        }

        let (suffix_array, reduced_text) =
            suffix_array.split_at_mut(suffix_array.len() - num_lms_substrings);
        let free_space = free_space + n - num_lms_substrings * 2;

        sais(reduced_text, suffix_array, name_count as usize, free_space);

        let mut j = num_lms_substrings;
        let mut c1 = text[n - 1];
        let mut s_type = false;

        for i in (0..n - 1).rev() {
            let c0 = text[i];
            if c0 < c1 {
                s_type = true;
            } else if c0 > c1 && s_type {
                reduced_text[j - 1] = I::from(i).unwrap() + I::ONE;
                j -= 1;
                s_type = false;
            }
            c1 = c0;
        }

        for i in 0..num_lms_substrings {
            suffix_array[i] = reduced_text[suffix_array[i].as_usize()];
        }
    }

    get_buckets(text, buckets, alphabet_size, true);
    buckets[separator.char_to_usize()] -= I::ONE;
    let separator_pos = buckets[separator.char_to_usize()].as_usize();
    suffix_array[num_lms_substrings..].fill(I::ZERO);
    for i in (0..num_lms_substrings).rev() {
        let j = suffix_array[i];
        let c = text[j.as_usize()].char_to_usize();
        suffix_array[i] = I::ZERO;
        buckets[c] -= I::ONE;
        suffix_array[buckets[c].as_usize()] = j;
    }

    suffix_array[separator_pos] = I::from(n - 1).unwrap();
    induce_sa_generalized(text, suffix_array, buckets, separator, alphabet_size, false);
}

pub fn gsais_lcp<C: Character, I: Index>(
    text: &[C],
    mut suffix_array: &mut [I],
    lcp: &mut [I],
    separator: C,
    alphabet_size: usize,
    free_space: usize,
) {
    let n = text.len();
    let buckets;
    if alphabet_size <= free_space {
        (suffix_array, buckets) = suffix_array.split_at_mut(suffix_array.len() - alphabet_size);
    } else {
        buckets = unsafe {
            let ptr = vec![I::ZERO; alphabet_size].as_mut_ptr();
            std::slice::from_raw_parts_mut(ptr, alphabet_size)
        }
    }

    let num_lms_substrings =
        sort_lms_substrings_generalized(text, suffix_array, separator, buckets, alphabet_size);
    let name_count =
        name_lms_substrings_generalized(text, suffix_array, num_lms_substrings, lcp, separator);

    let mut j = suffix_array.len() - 1;
    for i in (num_lms_substrings..n).rev() {
        if suffix_array[i] != I::ZERO {
            suffix_array[j] = suffix_array[i];
            j = j.saturating_sub(1);
        }
    }

    if (name_count as usize) < num_lms_substrings {
        let (suffix_array, reduced_text) =
            suffix_array.split_at_mut(suffix_array.len() - num_lms_substrings);

        let free_space = free_space + n - num_lms_substrings * 2;

        sais(reduced_text, suffix_array, name_count + 1, free_space);

        let mut j = num_lms_substrings;
        let mut c1 = text[n - 1];
        let mut s_type = false;

        for i in (0..n - 1).rev() {
            let c0 = text[i];
            if c0 < c1 {
                s_type = true;
            } else if c0 > c1 && s_type {
                reduced_text[j - 1] = I::from(i + 1).unwrap();
                j -= 1;
                s_type = false;
            }
            c1 = c0;
        }

        compute_lcp_phi_sparse(text, reduced_text, suffix_array, lcp, separator);

        for i in 0..num_lms_substrings {
            suffix_array[i] = reduced_text[suffix_array[i].as_usize()];
        }
    } else {
        for i in 1..num_lms_substrings {
            let s1 = text[suffix_array[i - 1].as_usize()..].iter();
            let s2 = text[suffix_array[i].as_usize()..].iter();

            let l = s1
                .zip(s2)
                .take_while(|(&a, &b)| a == b && !(a == separator && b == separator))
                .count();
            lcp[i] = I::from(l).unwrap();
        }
    }

    get_buckets(text, buckets, alphabet_size, true);
    buckets[separator.char_to_usize()] -= I::ONE;
    let separator_pos = buckets[separator.char_to_usize()].as_usize();

    suffix_array[num_lms_substrings..].fill(I::max_value());
    lcp[num_lms_substrings..].fill(I::ZERO);

    for i in (0..num_lms_substrings).rev() {
        let j = suffix_array[i];
        let l = lcp[i];
        let c = text[j.as_usize()].char_to_usize();

        suffix_array[i] = I::max_value();
        lcp[i] = I::ZERO;
        buckets[c] -= I::ONE;
        suffix_array[buckets[c].as_usize()] = j;
        lcp[buckets[c].as_usize()] = l;
    }

    suffix_array[separator_pos] = I::from(n - 1).unwrap();
    induce_sa_and_lcp_generalized(text, suffix_array, lcp, buckets, separator, alphabet_size);
}

fn sais<C: Character, I: Index>(
    text: &[C],
    mut suffix_array: &mut [I],
    alphabet_size: usize,
    free_space: usize,
) {
    let buckets;
    if alphabet_size <= free_space {
        (suffix_array, buckets) = suffix_array.split_at_mut(suffix_array.len() - alphabet_size);
    } else {
        buckets = unsafe {
            let ptr = vec![I::ZERO; alphabet_size].as_mut_ptr();
            std::slice::from_raw_parts_mut(ptr, alphabet_size)
        }
    }

    get_buckets(text, buckets, alphabet_size, true);
    let n = text.len();
    suffix_array.fill(I::ZERO);

    // place the LMS suffixes
    let mut c1 = text[n - 1];
    let mut s_type = false;

    for i in (0..n - 1).rev() {
        let c0 = text[i];
        if c0 < c1 {
            s_type = true;
        } else if c0 > c1 && s_type {
            buckets[c1.char_to_usize()] -= I::ONE;
            suffix_array[buckets[c1.char_to_usize()].as_usize()] = I::from(i + 1).unwrap();
            s_type = false;
        }
        c1 = c0;
    }

    // induce S and L type suffixes
    induce_sa(text, suffix_array, buckets, alphabet_size);

    // Store the sorted LMS-substrings in the first m items of SA
    let mut m: usize = 0;
    for i in 0..n {
        let pos = suffix_array[i];
        let c0 = text[pos.as_usize()];

        if pos > I::ZERO && text[pos.as_usize() - 1] > c0 {
            let mut j = pos.as_usize() + 1;

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
        name_buffer.fill(I::ZERO);

        let mut j = n;
        let mut c1 = text[j - 1];
        let mut is_stype = false;

        for i in (0..n - 1).rev() {
            let c0 = text[i];
            if c0 < c1 {
                is_stype = true;
            } else if c0 > c1 && is_stype {
                name_buffer[(i + 1) >> 1] = I::from(j - i - 1).unwrap();
                is_stype = false;
                j = i + 1;
            }
            c1 = c0;
        }

        let mut prev_pos = n;
        let mut prev_len = 0;

        for i in 0..m {
            let curr_pos = sa[i].as_usize();
            let curr_len = name_buffer[curr_pos >> 1].as_usize();
            let mut diff = true;

            if prev_len == curr_len {
                diff = text[curr_pos..curr_pos + curr_len] != text[prev_pos..prev_pos + prev_len];
            }

            if diff {
                name += 1;
                prev_pos = curr_pos;
                prev_len = curr_len;
            }

            name_buffer[curr_pos >> 1] = I::from(name).unwrap();
        }
    }

    if (name as usize) < m {
        let mut j = suffix_array.len() - 1;

        for i in (m..n).rev() {
            if suffix_array[i] != I::ZERO {
                suffix_array[j] = suffix_array[i] - I::ONE;
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
                reduced_text[j - 1] = I::from(i + 1).unwrap();
                j -= 1;
                s_type = false;
            }
            c1 = c0;
        }

        for i in 0..m {
            suffix_array[i] = reduced_text[suffix_array[i].as_usize()];
        }
    }

    get_buckets(text, buckets, alphabet_size, true);
    suffix_array[m..].fill(I::ZERO);
    for i in (0..m).rev() {
        let j = suffix_array[i];
        let c = text[j.as_usize()].char_to_usize();
        suffix_array[i] = I::ZERO;
        buckets[c] -= I::ONE;
        suffix_array[buckets[c].as_usize()] = j;
    }

    induce_sa(text, suffix_array, buckets, alphabet_size);
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
    use super::gsais;
    use super::gsais_lcp;

    #[test]
    fn try_suffix_array() {
        // let text = "TTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGGGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCATCTGGTAGCGATGATTGA";
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

    #[test]
    fn test_gsais() {
        // let text = "banana$bananabananabana$banananananananan$";
        //LSLSLLSSLSLLSSLSLL
        let text = "TTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCG$GGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCAC$CGTCCTCTCTGCCCCC$GCCAAAATCACCAACCATCTGGTAGCGATGATTGA$";
        // let text = "TTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGGGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCATCTGGTAGCGATGATTGA$";
        let mut suffix_array = vec![0; text.len() + 1];
        suffix_array[0] = text.len() as i32;
        gsais(text.as_bytes(), &mut suffix_array[1..], b'$', 256, 0);
        // let result = &[
        //     18, 6, 12, 17, 5, 11, 9, 15, 3, 7, 13, 1, 10, 0, 16, 4, 8, 14, 2
        // ];
        let result = &[
            144, 41, 90, 107, 143, 142, 111, 59, 112, 20, 120, 28, 60, 113, 21, 64, 88, 9, 117,
            121, 29, 26, 61, 131, 70, 74, 80, 14, 50, 114, 11, 124, 135, 22, 65, 52, 138, 89, 106,
            110, 19, 119, 63, 87, 116, 25, 69, 73, 79, 13, 10, 123, 105, 109, 118, 86, 68, 122,
            104, 103, 102, 47, 94, 39, 48, 133, 2, 5, 91, 30, 95, 97, 99, 36, 55, 126, 40, 141, 58,
            27, 49, 134, 137, 18, 62, 24, 72, 78, 108, 85, 67, 101, 46, 38, 132, 4, 57, 17, 71, 77,
            84, 3, 76, 75, 6, 128, 81, 42, 7, 129, 92, 15, 82, 43, 31, 8, 130, 51, 115, 12, 93, 1,
            96, 98, 35, 54, 125, 140, 136, 23, 66, 100, 45, 37, 56, 16, 83, 127, 0, 34, 53, 139,
            44, 33, 32,
        ];

        assert_eq!(&suffix_array, result);
    }

    #[test]
    fn test_gsais_lcp() {
        // let text = "banana$anaba$anan$";
        // let text = "banana$bananabananabana$banananananananan$";
        let text = "TTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCG$GGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCAC$CGTCCTCTCTGCCCCC$GCCAAAATCACCAACCATCTGGTAGCGATGATTGA$";
        let mut suffix_array = vec![0i16; text.len() + 1];
        let mut lcp = vec![0; suffix_array.len()];
        suffix_array[0] = text.len() as i16;
        gsais_lcp(
            text.as_bytes(),
            &mut suffix_array[1..],
            &mut lcp[1..],
            b'$',
            256,
            0,
        );

        let result = &[
            0, 0, 0, 0, 0, 0, 1, 3, 3, 4, 2, 3, 2, 2, 3, 5, 1, 2, 2, 4, 2, 1, 2, 3, 2, 3, 3, 2, 1,
            2, 4, 3, 2, 3, 4, 2, 3, 0, 1, 1, 4, 3, 3, 2, 3, 2, 3, 4, 4, 3, 2, 4, 1, 2, 4, 3, 3, 3,
            2, 3, 4, 2, 2, 1, 2, 4, 2, 3, 2, 3, 1, 4, 2, 4, 3, 4, 0, 1, 2, 3, 2, 3, 3, 1, 4, 3, 4,
            5, 2, 4, 4, 3, 3, 2, 3, 3, 1, 2, 4, 6, 3, 3, 2, 3, 2, 4, 3, 3, 1, 3, 2, 2, 5, 2, 3, 0,
            2, 2, 1, 3, 2, 2, 2, 3, 5, 4, 5, 1, 3, 2, 3, 4, 4, 3, 2, 3, 4, 3, 1, 3, 5, 2, 3, 2, 3,
        ];
        assert_eq!(&lcp, result);

        let result = &[
            144, 41, 90, 107, 143, 142, 111, 59, 112, 20, 120, 28, 60, 113, 21, 64, 88, 9, 117,
            121, 29, 26, 61, 131, 70, 74, 80, 14, 50, 114, 11, 124, 135, 22, 65, 52, 138, 89, 106,
            110, 19, 119, 63, 87, 116, 25, 69, 73, 79, 13, 10, 123, 105, 109, 118, 86, 68, 122,
            104, 103, 102, 47, 94, 39, 48, 133, 2, 5, 91, 30, 95, 97, 99, 36, 55, 126, 40, 141, 58,
            27, 49, 134, 137, 18, 62, 24, 72, 78, 108, 85, 67, 101, 46, 38, 132, 4, 57, 17, 71, 77,
            84, 3, 76, 75, 6, 128, 81, 42, 7, 129, 92, 15, 82, 43, 31, 8, 130, 51, 115, 12, 93, 1,
            96, 98, 35, 54, 125, 140, 136, 23, 66, 100, 45, 37, 56, 16, 83, 127, 0, 34, 53, 139,
            44, 33, 32,
        ];

        assert_eq!(&suffix_array, result);
    }
}
