use bbwt::parse::{PFBwt, Parse};
use bbwt::scan::Scanner;

use lexopt::prelude::*;
use seq_io::fasta;

use std::path::PathBuf;

#[derive(Debug)]
struct Args {
    window: u32,
    modulo: u32,
    output_path: PathBuf,
    filenames: Vec<String>,
}

fn parse_args() -> Result<Args, lexopt::Error> {
    let mut window = 10;
    let mut modulo = 100;
    let mut output_path = PathBuf::new();
    let mut filenames = Vec::new();
    let mut parser = lexopt::Parser::from_env();

    while let Some(arg) = parser.next()? {
        match arg {
            Short('w') | Long("window") => {
                window = parser.value()?.parse()?;
            }
            Short('m') | Long("modulus") => {
                modulo = parser.value()?.parse()?;
            }
            Short('o') | Long("output") => {
                output_path = parser.value()?.parse()?;
            }
            Value(filename) => {
                filenames.push(filename.string()?);
            }
            _ => return Err(arg.unexpected()),
        }
    }

    Ok(Args {
        window,
        modulo,
        output_path,
        filenames,
    })
}

fn main() {
    let mut args = parse_args().expect("Error processing command line arguments");
    args.output_path = if args.output_path.is_absolute() {
        args.output_path
    } else {
        let (dirname, basename) = (args.output_path.parent(), args.output_path.file_name());
        let mut canonical_dirname = match dirname {
            Some(path) => path.canonicalize().unwrap(),
            None => std::env::current_dir().unwrap(),
        };
        canonical_dirname.push(basename.unwrap());

        canonical_dirname
    };

    {
        let mut scanner = Scanner::new(&args.output_path);

        if args.filenames.is_empty() {
            let stdin = std::io::stdin();
            let reader = fasta::Reader::new(stdin.lock());
            scanner.process_fasta(reader, args.window, args.modulo);
        } else {
            for filename in args.filenames {
                let reader = fasta::Reader::from_path(filename).expect("Unable to open FASTA file");
                scanner.process_fasta(reader, args.window, args.modulo);
            }
        }
        scanner.serialize();
    }

    {
        let mut parse = Parse::new(&args.output_path);
        parse.transform_text_into_bwt();
        parse.contruct_inverse_list();
    }

    let mut pfbwt = PFBwt::new(&args.output_path);
    pfbwt.contruct_bwt(args.window);
}
