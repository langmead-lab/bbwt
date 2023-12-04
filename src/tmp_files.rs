use std::io::{BufRead, Write};

struct TempFiles {
    last_file: String
}

enum Mode<'a> {
    Reader(&'a dyn BufRead),
    Writer(&'a dyn Write)
}

impl TempFiles {
    fn new(basename: String) -> Self {
        unimplemented!()
    }

    fn open_last_file(&mut self) {
    }
}
