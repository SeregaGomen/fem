extern crate json;
use std::env;
use fem::json::read_json; 

fn main() {
    let args: Vec<String> = env::args().collect();
   
    if args.len() < 2 {
        println!("Too few parameters!");
    } else {
        if let Err(e) = read_json(args[1].as_str()) {
            // println!("\n\x1b[93m{}\x1b[0m", e);
            println!("\n\x1b[91m{}\x1b[0m", e);
        }
    }
}

