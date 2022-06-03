extern crate json;

use fem::json::read_json; 

fn main() {
    use std::env;

    let args: Vec<String> = env::args().collect();
    
    if args.len() < 2 {
        println!("Too few parameters!");
        return;
    }
    match read_json(args[1].as_str()) {
         Ok(_) => println!("Done"),
         Err(e) => println!("Error: {}", e),
    };
}

