extern crate bio;
extern crate rand;
extern crate suffix_tree;

use suffix_tree::SuffixTree;
use bio::alphabets;
use std::io::{self, Write};
use std::time::{Instant};

fn main() {
    let mut input = String::new();
    print!("Enter suffix tree string: ");
    io::stdout().flush().unwrap();
    match io::stdin().read_line(&mut input) {
        Err(_) => return,
        Ok(_) => ()
    }

    let now = Instant::now();
    let tree = SuffixTree::new(alphabets::dna::alphabet(), &input.trim().to_string());
    let elapsed = now.elapsed();
    println!("Suffix tree created in {}.{} seconds.", elapsed.as_secs(), elapsed.subsec_nanos());

    let mut search = String::new();
    print!("Enter search string: ");
    io::stdout().flush().unwrap();
    match io::stdin().read_line(&mut search) {
        Err(_) => return,
        Ok(_) => ()
    }

    print!("The full string ");
    if tree.contains(&search.trim().to_string()) {
        print!("contains");
    } else {
        print!("does not contain");
    }
    println!(" the search string.");
}
