use rayon::prelude::*;
use ndarray::Array1;
use crate::error::Error;

pub trait SparseMatrix {
    fn add_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error>;
    fn set_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error>;
    fn get_value(&mut self, index1: usize, index2: usize) -> Result<f64, Error>;
    fn clear_row(&mut self, index: usize) -> Result<(), Error>;
    fn clear_col(&mut self, index: usize) -> Result<(), Error>;
}

pub struct MapSparseMatrix {
    nvtxs: usize,
    blksze: usize,
    map: Vec<Vec<usize>>,
    data: Vec<Vec<f64>>,
}

impl SparseMatrix for MapSparseMatrix {
    fn add_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error> {
        let pos = self.find(index1, index2)?;
        //println!("{} - {}", pos.0, pos.1);
        self.data[pos.0][pos.1] += value;
        Ok(())
    }
    fn set_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error> {
        let pos = self.find(index1, index2)?;
        self.data[pos.0][pos.1] = value;
        Ok(())
    }
    fn get_value(&mut self, index1: usize, index2: usize) -> Result<f64, Error> {
        let pos = self.find(index1, index2)?;
        Ok(self.data[pos.0][pos.1])
    }
    fn clear_row(&mut self, index: usize) -> Result<(), Error> {
        if index >= self.nvtxs * self.blksze {
            return Err(Error::InvalidIndex);    
        }
        for i in 0..self.data[index].len() {
            self.data[index][i] = 0.;
        }
        Ok(())
    }
    fn clear_col(&mut self, index: usize) -> Result<(), Error> {
        if index >= self.nvtxs * self.blksze {
            return Err(Error::InvalidIndex);    
        }
        for i in 0..self.nvtxs * self.blksze {
            match self.find(i, index) {
                Err(_) => continue,
                Ok(pos) => {
                    self.data[i][pos.1] = 0.;
                }
                    
            }
        }
        Ok(())
    }

}


#[allow(dead_code)]
impl MapSparseMatrix {
    pub fn new(nvtxs: usize, blksze: usize, map: &Vec<Vec<usize>>) -> Self {
        let mut data: Vec<Vec<f64>> = Vec::new();
        for i in 0..nvtxs * blksze {
            data.push(vec![0.0; map[i / blksze].len() * blksze]);
        }
        Self{nvtxs, blksze, map: map.clone(), data}
    }
    fn find(&self, index1: usize, index2: usize) -> Result<(usize, usize), Error> {
        let row: usize = index1 / self.blksze;
        if index1 >= self.nvtxs * self.blksze || index2 >= self.nvtxs * self.blksze  {
            return Err(Error::InvalidIndex);
        } 
        for i in 0..self.map[row].len() {
            if self.map[row][i] == index2 / self.blksze  {
                return Ok((index1, i * self.blksze + index2 % self.blksze));
            }
        }
        Err(Error::InvalidIndex)
    }
    // pub fn dot(&self, rhs: &Array1<f64>) -> Array1<f64> {
    //     let mut x: Array1<f64> = Array1::zeros(self.size * self.freedom);
    //     for i in 0..self.rows() {
    //         for j in 0..self.data[i].len() {
    //             let col = self.map[i / self.freedom][j / self.freedom] * self.freedom + j % self.freedom;
    //             x[i] += self.data[i][j] * rhs[col];
    //         }
    //     }
    //     x
    // }
    pub fn dot(&self, rhs: &Array1<f64>) -> Array1<f64> {
        let mut rows = (0..self.nvtxs * self.blksze)
            .into_par_iter()
            .map(move |i| {
                (i, (0..self.data[i].len())
                .map(|j| {
                    let col = self.map[i / self.blksze][j / self.blksze] * self.blksze + j % self.blksze;
                    &self.data[i][j] * &rhs[col]
                }).sum::<f64>())
            })
            .collect::<Vec<(usize, f64)>>();
        rows.par_sort_by(|left, right| left.0.cmp(&right.0));
        rows.into_iter().map(|(_, row)| row).collect::<Array1<f64>>()
    }
}
