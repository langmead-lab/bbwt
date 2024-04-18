use crate::int::Index;
use crate::int::Character;

#[derive(Copy, Clone, Debug)]
struct LCPEntry<I: Index> {
    index: usize,
    lcp: I,
}

pub(crate) struct StackRMQ<I: Index> {
    stack: Vec<LCPEntry<I>>,
    last_occ: Vec<usize>,
    sorted_occ: Vec<usize>,
    s_type: bool,
    cmp: fn(&usize, &usize) -> std::cmp::Ordering,
}

impl<I: Index> StackRMQ<I> {
    pub(crate) fn new(n: usize, s_type: bool) -> Self {
        let stack = Vec::with_capacity(2 * n);
        let empty_occ = if s_type { usize::MIN } else { usize::MAX };
        let last_occ = vec![empty_occ; n];
        let sorted_occ = vec![0; n];

        StackRMQ {
            stack,
            last_occ,
            sorted_occ,
            s_type,
            cmp: if !s_type {
                |a, b| a.cmp(b)
            } else {
                |a, b| b.cmp(a)
            },
        }
    }

    pub(crate) fn reset_for_type(&mut self, s_type: bool) {
        let empty_occ = if s_type { usize::MIN } else { usize::MAX };
        self.s_type = s_type;
        self.last_occ.fill(empty_occ);
        self.stack.truncate(0);
        self.cmp = if !s_type {
            |a, b| a.cmp(b)
        } else {
            |a, b| b.cmp(a)
        };
    }

    pub(crate) fn push(&mut self, index: usize, lcp: I) {
        self.remove_stale_lcps(lcp);
        self.stack.push(LCPEntry { index, lcp });
        dbg!(&self.stack);
    }

    pub(crate) fn get_min<C: Character>(&mut self, character: C, index: usize) -> I {
        let mut min_lcp = I::ZERO;
        let empty_occ = if self.s_type { usize::MIN } else { usize::MAX };
        let last = self.get_and_update_last_occ(character, index);

        if last == empty_occ {
            return min_lcp;
        }

        if let Some(entry) = self.stack.iter().find(|e| {
            if self.s_type {
                e.index <= last
            } else {
                e.index > last
            }
        }) {
            min_lcp = entry.lcp + I::ONE;
        }

        if self.stack.len() > self.last_occ.len() {
            self.shrink_stack();
        }

        min_lcp
    }

    #[inline]
    fn get_and_update_last_occ<C: Character>(&mut self, character: C, pos: usize) -> usize {
        let last = self.last_occ[character.char_to_usize()];
        self.last_occ[character.char_to_usize()] = pos;

        last
    }

    #[inline]
    fn remove_stale_lcps(&mut self, lcp: I) {
        let new_size =
            self.stack.len() - self.stack.iter().rev().take_while(|e| e.lcp >= lcp).count();

        self.stack.truncate(new_size);
    }

    #[inline]
    fn shrink_stack(&mut self) {
        let mut cur = 0;
        let mut end = 0;
        let empty_occ = if self.s_type { usize::MAX } else { usize::MIN };
        self.sorted_occ.extend_from_slice(&self.last_occ);
        self.sorted_occ.sort_by(self.cmp);

        for &occ in &self.sorted_occ {
            if occ == empty_occ {
                continue;
            }
            if self.stack[end].index < occ {
                cur += self.stack[cur..]
                    .iter()
                    .take_while(|&e| e.index < occ)
                    .count()
                    + 1;
            }
            if cur < self.stack.len() {
                self.stack[end] = self.stack[cur];
                end += 1;
                cur += 1;
            }
        }

        self.stack.truncate(end);
        self.sorted_occ.truncate(0);
    }
}

pub(crate) fn compute_lcp_phi_sparse<C: Character, I: Index>(
    text: &[C],
    reduced_text: &[I],
    suffix_array: &mut [I],
    lcp: &mut [I],
    separator: C,
) {
    let n1 = reduced_text.len();
    let (lcp, plcp) = lcp.split_at_mut(n1);

    plcp[suffix_array[0].as_usize()] = I::ZERO;
    for i in 1..n1 {
        plcp[suffix_array[i].as_usize()] = lcp[i];
    }

    lcp[suffix_array[0].as_usize()] = I::ZERO;
    for i in 1..n1 {
        lcp[suffix_array[i].as_usize()] = suffix_array[i - 1];
    }

    let mut l = I::ZERO;
    for i in 0..n1 - 1 {
        if text[reduced_text[i].as_usize()] == separator {
            continue;
        }

        l = plcp[i].max(l);
        let s1 = reduced_text[i] + l;
        let s2 = reduced_text[lcp[i].as_usize()] + l;

        l += I::from(text[s1.as_usize()..]
            .iter()
            .copied()
            .zip(text[s2.as_usize()..].iter().copied())
            .take_while(|(a, b)| a == b && !(*a == separator && *b == separator))
            .count()).unwrap();
        plcp[i] = I::from(l).unwrap();

        l -= if lcp[i].as_usize() == n1 - 1 {
            reduced_text[i + 1] - reduced_text[i]
        } else {
            (reduced_text[i + 1] - reduced_text[i])
                .max(reduced_text[lcp[i].as_usize() + 1] - reduced_text[lcp[i].as_usize()])
        }
    }

    lcp[0] = I::ZERO;
    for i in 1..n1 {
        lcp[i] = plcp[suffix_array[i].as_usize()];
    }
}

#[cfg(test)]
mod tests {
    use super::compute_lcp_phi_sparse;

    #[test]
    fn try_compute_lcp() {
        let text = "banaananaanana$";
        let mut sa = vec![5, 3, 1, 4, 2, 0, 0, 0, 0, 1, 3, 6, 8, 11, 14];
        let len = sa.len();
        let mut lcp = vec![0, 0, 4, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let (suffix_array, reduced_text) = sa.split_at_mut(len - 6);

        compute_lcp_phi_sparse(text.as_bytes(), reduced_text, suffix_array, &mut lcp, b'$');

        assert_eq!(lcp, &[0, 0, 6, 1, 3, 8, 8, 6, 3, 0, 1, 0, 0, 0, 0]);
    }
}
