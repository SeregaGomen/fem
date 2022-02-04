use rayon::prelude::*;
use ndarray::{Array1, prelude::*};
use crate::error::Error;
use super::util;
use super::msg::Messenger;


pub trait SparseMatrix {
    fn add_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error>;
    fn set_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error>;
    fn get_value(&self, index1: usize, index2: usize) -> Result<f64, Error>;
    fn solve(&mut self, rhs: &Array1<f64>, eps: f64) -> Result<Array1<f64>, Error>;
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
    size: usize,
    diag: Vec<f64>,
    env: Vec<f64>,
    xenv: Vec<usize>,
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
    fn get_value(&self, index1: usize, index2: usize) -> Result<f64, Error> {
        let pos = self.find(index1, index2)?;
        Ok(self.data[pos.0][pos.1])
    }
    fn solve(&mut self, rhs: &Array1<f64>, eps: f64) -> Result<Array1<f64>, Error> {
        self.lzh_solve(rhs, eps)
    }
}


impl SparseMatrix for EnvSparseMatrix {
    fn add_value(&mut self, i: usize, j: usize, value: f64) -> Result<(), Error> {
        if i >= j {
            if i == j {
                self.diag[i] += value;
            }
            else
            {
                if self.xenv[i + 1] - i + j >= self.xenv[i] {
                    self.env[self.xenv[i + 1] - i + j] += value;
                }
            }
        }
        Ok(())
    }
    fn set_value(&mut self, i: usize, j: usize, value: f64) -> Result<(), Error> {
        if i >= j {
            if i == j {
                self.diag[i] = value;
            }
            else
            {
                if self.xenv[i + 1] - i + j >= self.xenv[i] {
                    self.env[self.xenv[i + 1] - i + j] = value;
                }
            }
        }
        Ok(())
    }
    fn get_value(&self, i: usize, j: usize) -> Result<f64, Error> {
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
    fn solve(&mut self, rhs: &Array1<f64>, _eps: f64) -> Result<Array1<f64>, Error> {
        // Факторизация матрицы A = L * L(t)
        if self.esfct() == false {
            return Err(Error::SingularMatrix);
        }
        // Вычисление y: (Ly = b)
        let y = self.elslv(self.size, &self.diag, &self.env, &self.xenv, &rhs.to_vec());
        // Вычисление x: (L(t)x = y)
        let x = self.euslv(self.size, &self.diag, &self.env, &self.xenv, &y);
        Ok(Array::from_vec(x))
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
    // fn dot(&self, rhs: &Array1<f64>) -> Array1<f64> {
    //     let mut x: Array1<f64> = Array1::zeros(self.size * self.freedom);
    //     for i in 0..self.rows() {
    //         for j in 0..self.data[i].len() {
    //             let col = self.map[i / self.freedom][j / self.freedom] * self.freedom + j % self.freedom;
    //             x[i] += self.data[i][j] * rhs[col];
    //         }
    //     }
    //     x
    // }
    fn dot(&self, rhs: &Array1<f64>) -> Array1<f64> {
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
    // Метод Ланцоша
    fn lzh_solve(&self, rhs: &Array1<f64>, eps: f64) -> Result<Array1<f64>, Error> {
        let max_iter = 2 * self.nvtxs * self.blksze;
        let mut x: Array1<f64> = rhs.clone();
        let mut r = self.dot(&x) - &x;
        let mut s = r.clone();
        let mut norm = util::scalar_product(&r, &r);
        let mut is_ok = false;
        
        // println!("{:?}", x);
        let mut msg = Messenger::new("Solution of the system of equations", 1, max_iter as i64, 1);
        for _i in 0..max_iter {
            msg.add_progress();
            let r1 = self.dot(&s);
            let err = util::scalar_product(&r1, &s);
            // if i % 100 == 0 {
            //     println!("Error: {}", err);
            // }
            if err.abs() < eps {
                // println!("Error: {}", err);
                msg.stop();
                is_ok = true;
                break;
            }
            let mut alpha = norm / err;
            r = &r - &r1 * alpha;
            x = &x - &s * alpha;
            let norm1 = util::scalar_product(&r, &r);
            alpha = norm1 / norm;
            s = &s * alpha + &r;
            norm = norm1;
        }
        if !is_ok {
            return Err(Error::SingularMatrix);
        }
        Ok(x)
    }
}

#[allow(dead_code)]
impl EnvSparseMatrix {
    pub fn new(nvtxs: usize, blksze: usize, map: &Vec<Vec<usize>>) -> Self {
        let mut size: usize = 0;
        let mut count: usize = 0;
        let diag = vec![0.; nvtxs * blksze];
        let env: Vec<f64>;
        let mut xenv = vec![0; nvtxs * blksze + 1];
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
        env = vec![0.; size];
        size = xenv.len() - 1;
        xenv[size] = env.len();
        Self{size: nvtxs * blksze, diag, env, xenv}
    }
    fn esfct(&mut self) -> bool {
        if self.diag[0] < 0. {
            return false;
        }
        self.diag[0] = self.diag[0].sqrt();
        if self.size == 1 {
            return true;
        }
        let mut msg = Messenger::new("Factorization of the system of equations", 1, self.size as i64, 1);
        for i in 1..self.size {
            msg.add_progress();
            let ixenv = self.xenv[i];
            let iband = self.xenv[i + 1] - ixenv;
            let mut temp = self.diag[i];
            if iband != 0 {
                let ifirst = i - iband;
                // Вычислить строку i треугольного множителя
                ////////////////
                let res = self.elslv(iband, &self.diag[ifirst..], &self.env, &self.xenv[ifirst..], &self.env[ixenv..]);
                // env[ifirst..] = res;   
                for j in 0..res.len() {
                    self.env[ixenv + j] = res[j];
                }         
                ////////////////
    
    
                let jstop = self.xenv[i + 1];
                for j in ixenv..jstop {
                    let s = self.env[j]; 
                    temp -= s * s;
                }
            }
            if temp <= 0. {
                // Матрица не положительно определенная
                return false;
            }
            self.diag[i] = temp.sqrt();
        }
        msg.stop();
        true
    }
    // Решение нижней треугольной системы L*X = RHS
    fn elslv(&self, neqns: usize, diag: &[f64], env: &[f64], xenv: &[usize], rhs: &[f64]) -> Vec<f64> {
        let mut res = vec![0.; neqns];
        let mut ifirst: usize = 0;
        let mut last: usize = 0;
        // Найти номер первого ненулевого элемента в правой части и поместить его в ifirst
        for i in 0..neqns {
            if rhs[i] != 0. {
                ifirst = i;
                break;
            }
        }
        // last содержит номер последней вычисленной ненулевой компоненты решения
        for i in ifirst..neqns {
            let mut iband = xenv[i + 1] - xenv[i];
            if iband > i {
                iband = i;
            } 
            let mut s = rhs[i];
            let mut l = i - iband;
            // Строка оболочки пуста или все соответствующие компоненты решения - нули
            if !(iband == 0 || last < l) {
                let kstrt = xenv[i + 1] - iband;
                let kstop = xenv[i + 1];
                for k in kstrt..kstop {
                    s -= env[k] * res[l];
                    l += 1;
                }
            }
            if s != 0. {
                res[i] = s / diag[i];
                last = i;
            }
        }
        res
    }

    // Решение верхней треугольной системы профильным методом
    fn euslv(&self, neqns: usize, diag: &[f64], env: &[f64], xenv: &[usize], rhs: &[f64]) -> Vec<f64> {
        let mut res: Vec<f64> = Vec::new();
        let mut msg = Messenger::new("Solution of the system of equations", 1, neqns as i64, 1);
        res.extend_from_slice(rhs);
        for i in (0..neqns).rev() {
            msg.add_progress();
            if res[i] != 0. {
                let s = res[i] / diag[i];
                res[i] = s;
                let mut iband = xenv[i + 1] - xenv[i];
                if iband > i {
                    iband -= 1;
                }
                if iband == 0 {
                    continue;
                }
                let kstrt = i - iband;
                let kstop = i;
                let mut l = xenv[i + 1] - iband;
                for k in kstrt..kstop {
                    res[k] -= s * env[l];
                    l += 1; 
                }
            }
        }
        res
    }
}