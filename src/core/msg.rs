use std::io::{self, Write};
use std::time::Instant;

pub struct Messenger {
    msg: String,
    start: i64,
    stop: i64,
    current: i64,
    step: i64,
    old: i64,
    time: Instant,
}

impl Messenger {
    pub fn new(msg: String, start: i64, stop: i64, step: i64) -> Self {
        print!("\r{}...{}%", msg, 0);
        io::stdout().flush().unwrap();
        Self { msg, start, stop, step, current: 0, old: 0, time: Instant::now(), }
    }
    pub fn add_progress(&mut self) {
        self.current += 1;
        let persent: i64 = if self.stop - self.start != 0 { ((100.0 * self.current as f64) / ((self.stop - self.start) as f64)) as i64 } else { 100 };
        if self.current == self.stop {
            //let end = Instant::now();
            print!("\r{}...100%\n", self.msg);
            println!("Done in: {:.2?}", self.time.elapsed());
            return;
        }
        if persent == self.old {
            return;
        }
        if persent % self.step == 0 {
            print!("\r{}...{}%", self.msg, persent);
            io::stdout().flush().unwrap();
        }
        self.old = persent;
    }
    pub fn stop(&mut self) {
        print!("\r{}...100%\n", self.msg);
        println!("Done in: {:.2?}", self.time.elapsed());
    }
}