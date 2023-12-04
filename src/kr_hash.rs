pub(crate) struct KarpRabinHash {
    window_size: u32,
    asize: u64,
    prime: u64,
    hash: u64,
    char_count: u64,
    asize_pot: u64,
    window: Vec<u8>,
}

impl KarpRabinHash {
    pub(crate) fn new(window_size: u32) -> Self {
        let prime: u64 = 1999999973;
        let asize = 256;
        let mut asize_pot = 1u64;
        let window: Vec<u8> = vec![0; window_size as usize];

        for _ in 1..window_size {
            asize_pot = (asize_pot * asize as u64) % prime;
        }

        Self {
            window_size,
            asize,
            prime,
            hash: 0,
            char_count: 0,
            asize_pot,
            window,
        }
    }

    pub(crate) fn add_character(&mut self, char: u8) -> u64 {
        let index = (self.char_count % self.window_size as u64) as usize;
        self.char_count += 1;
        self.hash += self.prime - (self.window[index] as u64 * self.asize_pot) % self.prime;
        self.hash = (self.asize * self.hash + char as u64) % self.prime;
        self.window[index] = char;

        self.hash
    }

    pub(crate) fn kr_hash(b: impl AsRef<[u8]>) -> u64 {
        let mut hash: u64 = 0;
        let prime: u64 = 27162335252586509;

        for c in b.as_ref() {
            hash = (256 * hash + *c as u64) % prime;
        }

        hash
    }
}
