use std::io::Write;
use ndarray::Array1;
use crate::error::Error;

pub struct SparseMatrix {
    size: usize,
    freedom: usize,
    index: Vec<Vec<usize>>,
    data: Vec<Vec<f64>>,
}

#[allow(dead_code)]
impl SparseMatrix {
    pub fn new(size: usize, freedom: usize, ind: &Vec<Vec<usize>>) -> Self {
        let mut data: Vec<Vec<f64>> = Vec::new();
        for i in 0..size * freedom {
            data.push(vec![0.0; ind[i / freedom].len() * freedom]);
        }
        Self{size, freedom, index: ind.clone(), data}
    }
    pub fn rows(&self) -> usize {
        self.size * self.freedom
    }
    pub fn cols(&self) -> usize {
        self.size * self.freedom
    }
    fn find(&self, index1: usize, index2: usize) -> Result<(usize, usize), Error> {
        let row: usize = index1 / self.freedom;
        if index1 >= self.rows() || index2 >= self.cols()  {
            return Err(Error::InvalidIndex);
        } 
        for i in 0..self.index[row].len() {
            if self.index[row][i] == index2 / self.freedom  {
                return Ok((index1, i * self.freedom + index2 % self.freedom));
            }
        }
        Err(Error::InvalidIndex)
    }
    pub fn add_element(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error> {
        let pos = self.find(index1, index2)?;
        //println!("{} - {}", pos.0, pos.1);
        self.data[pos.0][pos.1] += value;
        Ok(())
    }
    pub fn set_element(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error> {
        let pos = self.find(index1, index2)?;
        self.data[pos.0][pos.1] = value;
        Ok(())
    }
    pub fn get_element(&mut self, index1: usize, index2: usize) -> Result<f64, Error> {
        let pos = self.find(index1, index2)?;
        Ok(self.data[pos.0][pos.1])
    }
    pub fn print(&mut self, file_name: &str, eps: f64) -> Result<(), Error> {
        let size = self.size * self.freedom;
        let mut file = match std::fs::File::create(file_name) {
            Err(_) => return Err(Error::OpenFile),
            Ok(file) => file,
        };
        file.write(size.to_string().as_bytes()).expect("write failed");
        file.write("\n".as_bytes()).expect("write failed");
        for i in 0..size {
            for j in 0..size {
                if self.get_element(i, j).is_err() {
                    file.write("0 ".as_bytes()).expect("write failed");
                }
                else {
                    if self.get_element(i, j)?.abs() > eps {
                        file.write("x ".as_bytes()).expect("write failed");
                    }
                    else {
                        file.write("0 ".as_bytes()).expect("write failed");
                    }

                }
            }
            file.write("\n".as_bytes()).expect("write failed");
        }
        Ok(())
    } 
    pub fn dot(&self, rhs: &Array1<f64>) -> Array1<f64> {
        let mut x: Array1<f64> = Array1::zeros(self.size * self.freedom);
        for i in 0..self.rows() {
            for j in 0..self.data[i].len() {
                let col = self.index[i / self.freedom][j / self.freedom] * self.freedom + j % self.freedom;
                x[i] += self.data[i][j] * rhs[col];
            }
        }
        x
    }
}
