use rayon::prelude::*;
use ndarray::{Array1, prelude::*};
use crate::error::Error;

pub trait SparseMatrix {
    fn add_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error>;
    fn set_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error>;
    fn get_value(&mut self, index1: usize, index2: usize) -> Result<f64, Error>;
    fn clear_row(&mut self, index: usize) -> Result<(), Error>;
    fn clear_col(&mut self, index: usize) -> Result<(), Error>;
}

// Разреженная матрица, хранящая только ненулевые элементы
pub struct MapSparseMatrix {
    nvtxs: usize,
    blksze: usize,
    map: Vec<Vec<usize>>,
    data: Vec<Vec<f64>>,
}

// Профильная разреженная матрица
pub struct EnvSparseMatrix {
    nvtxs: usize,
    blksze: usize,
    diag: Array1<f64>,
    env: Array1<f64>,
    xenv: Array1<usize>,
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
                Ok(pos) => self.data[i][pos.1] = 0.,
            }
        }
        Ok(())
    }
}


impl SparseMatrix for EnvSparseMatrix {
    fn add_value(&mut self, i: usize, j: usize, value: f64) -> Result<(), Error> {
        if i < j {
            return Err(Error::InvalidIndex);
        }
        if i == j {
            self.diag[i] += value;
        }
        else {
            if self.xenv[i + 1] - i + j >= self.xenv[i] {
                self.env[self.xenv[i + 1] - i + j] += value;            
            }
            else {
                return Err(Error::InvalidIndex);
            }
        }
        Ok(())
    }
    fn set_value(&mut self, i: usize, j: usize, value: f64) -> Result<(), Error> {
        if i < j {
            return Err(Error::InvalidIndex);
        }
        if i == j {
            self.diag[i] = value;
        }
        else {
            if self.xenv[i + 1] - i + j >= self.xenv[i] {
                self.env[self.xenv[i + 1] - i + j] = value;            
            }
            else {
                return Err(Error::InvalidIndex);
            }
        }
        Ok(())
    }
    fn get_value(&mut self, i: usize, j: usize) -> Result<f64, Error> {
        if i >= j {
            if i == j {
                return Ok(self.diag[i]);
            }
            else {
                if self.xenv[i + 1] - i + j >= self.xenv[i] {
                    return Ok(self.env[self.xenv[i + 1] - i + j]);
                }
            }
        }
        Err(Error::InvalidIndex)
    }
    fn clear_row(&mut self, index: usize) -> Result<(), Error> {
        if index >= self.nvtxs * self.blksze {
            return Err(Error::InvalidIndex);    
        }
        self.diag[index] = 0.;
        for i in self.xenv[index]..self.xenv[index + 1] - self.xenv[index] {
            self.env[i] = 0.;
        }
        Ok(())
    }
    fn clear_col(&mut self, index: usize) -> Result<(), Error> {
        if index >= self.nvtxs * self.blksze {
            return Err(Error::InvalidIndex);    
        }
        for i in 0..self.nvtxs * self.blksze {
            self.set_value(i, index, 0.)?;                    
        }
        Ok(())
    }
}

//#[allow(dead_code)]
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

impl EnvSparseMatrix {
    pub fn new(nvtxs: usize, blksze: usize, map: &Vec<Vec<usize>>) -> Self {
        let mut size: usize = 0;
        let mut count: usize = 0;
        let diag = Array1::<f64>::zeros(nvtxs * blksze);
        let env: Array1<f64>;
        let mut xenv = Array1::<usize>::zeros(nvtxs * blksze + 1);
        for i in 0..nvtxs {
            // Длина профиля для i-го узла
            let len = blksze * blksze * (i - map[i][0]) + (blksze - 1) * blksze / 2;
            if len == 0 {
                continue;
            }
            size += len;
            for k in 0..blksze {
                xenv[i * blksze + k] = count;
                count += blksze * (i - map[i][0]) + k;
            }
        }
        env = Array1::zeros(size);
        size = xenv.len() - 1;
        xenv[size] = env.len();
        Self{nvtxs, blksze, diag, env, xenv}
    }
}