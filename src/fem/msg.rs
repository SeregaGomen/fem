use std::io::{self, Write};
use std::time::Instant;
use std::sync::mpsc::{Sender, Receiver, channel, TryRecvError};

pub struct Messenger<'a> {
    msg: &'a str,
    start: i64,
    stop: i64,
    current: i64,
    step: i64,
    old: i64,
    time: Instant,
    tx: Sender<bool>,
}

impl<'a> Messenger<'a> {
    pub fn new(msg: &'static str, start: i64, stop: i64, step: i64) -> Self {
        let (tx, rx): (Sender<bool>, Receiver<bool>) = channel();
        if stop == 0 {
            // Перманентный процесс
            print!("\r{}...", msg);
            io::stdout().flush().unwrap();
            std::thread::spawn(move || {
                let mut i = 0;
                let chr = ['|', '/', '-', '\\'];
                loop {
                    print!("\r{}...{}", msg, chr[i % 4]);    
                    i += 1;
                    match rx.try_recv() {
                        Ok(_) | Err(TryRecvError::Disconnected) => break,
                        Err(TryRecvError::Empty) => {}
                    }
                    std::thread::sleep(std::time::Duration::from_millis(10));
                }
            });            
        } else {
            print!("\r{}...{}%", msg, 0);
            io::stdout().flush().unwrap();
        }
        Self { msg, start, stop, step, current: 0, old: 0, time: Instant::now(), tx }
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
        if self.stop == 0 {
            self.tx.send(true).unwrap();
            std::thread::sleep(std::time::Duration::from_millis(200));
        }
        print!("\r{}...100%\n", self.msg);
        println!("Done in: {:.2?}", self.time.elapsed());
    }
}