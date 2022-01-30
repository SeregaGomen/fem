use ndarray::{Array1};
use crate::error::Error;
use super::sparse::{SparseMatrix, MapSparseMatrix};
use super::mesh::Mesh;
use super::msg::Messenger;
use super::util;

pub trait Solver {
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error>;
    fn set_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error>;
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), Error>; 
    fn set_vector_value(&mut self, index: usize, val: f64) -> Result<(), Error>;
    fn set_result_value(&mut self, index: usize, val: f64) -> Result<(), Error>;
}

pub struct LzhSolver {
    nvtxs: usize,
    blksze: usize,
    a: MapSparseMatrix,
    b: Array1<f64>,
}

#[allow(dead_code)]
impl LzhSolver {
    pub fn new(mesh: &Mesh) -> Self {
        Self { 
            nvtxs: mesh.num_vertex,
            blksze: mesh.freedom,
            a: MapSparseMatrix::new(mesh.num_vertex, mesh.freedom, &mesh.mesh_map), 
            b: Array1::zeros(mesh.num_vertex * mesh.freedom), 
        }
    }
    // Метод сопряженных градиентов
    pub fn cg_solve(&self, eps: f64) -> Result<Array1<f64>, Error> {
        let max_iter = 100 * self.nvtxs * self.blksze;
        let mut x = self.b.clone();
        let mut r = &x - self.a.dot(&x);
        let mut z = r.clone();
        let mut is_ok = false;
        
        // self.a.print();
        // println!();
        let mut msg = Messenger::new("Solution of the system of equations", 1, 2 * max_iter as i64, 1);
        for _i in 0..max_iter {
            msg.add_progress();
            let p = self.a.dot(&z);
            let alpha = util::scalar_product(&r, &r) / util::scalar_product(&p, &z);
            x = &x + alpha * &z;  
            let r1 = &r - alpha * &p;
            let beta = util::scalar_product(&r1, &r1) / util::scalar_product(&r, &r);
            r = r1;
            if beta.sqrt() < eps {
                // println!("error: {}", beta.sqrt());
                msg.stop();
                is_ok = true;
                break;
            }
            z = &r + beta * &z;
            // if i % 100 == 0 {
            //     println!("error: {}", beta.sqrt());
            // }
        }
        if !is_ok {
            return Err(Error::SingularMatrix);
        }
        Ok(x)
    }
    // Метод Ланцоша
    pub fn lzh_solve(&self, eps: f64) -> Result<Array1<f64>, Error> {
        let max_iter = 2 * self.nvtxs * self.blksze;
        let mut x = self.b.clone();
        let mut r = self.a.dot(&x) - &x;
        let mut s = r.clone();
        let mut norm = util::scalar_product(&r, &r);
        let mut is_ok = false;
        
        // println!("{:?}", x);
        let mut msg = Messenger::new("Solution of the system of equations", 1, max_iter as i64, 1);
        for _i in 0..max_iter {
            msg.add_progress();
            let r1 = self.a.dot(&s);
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
    // Метод Холеского
    // Нижняя треугольная матрица L представлена в виде вектора длиной (n + 1) * n / 2
    // L{ij} -> l{(i + 1) * i / 2 + j}, i <= j
    pub fn cho_solve(&mut self, eps: f64) -> Result<Array1<f64>, Error> {
        let n = self.nvtxs * self.blksze;
        let mut l: Array1<f64> = Array1::zeros((n + 1) * n / 2);
        let mut msg = Messenger::new("Factorization of matrix", 1, n as i64, 1);
        // Факторизация матрицы A
        for i in 0..n {
            msg.add_progress();
            l[(i + 1) * i / 2 + i] = self.a.get_value(i, i)?;
            for k in 0..i {
                if l[(i + 1) * i / 2 + k] != 0. {
                    l[(i + 1) * i / 2 + i] -= l[(i + 1) * i / 2 + k] * l[(i + 1) * i / 2 + k];
                }
            }
            l[(i + 1) * i / 2 + i] = l[(i + 1) * i / 2 + i].sqrt();
            if l[(i + 1) * i / 2 + i].abs() < eps {
                return Err(Error::SingularMatrix);
            }
            for j in i + 1..n {
                l[(j + 1) * j / 2 + i] = match self.a.get_value(j, i) {
                    Ok(v) => v,
                    Err(_) => 0.,
                };
                for k in 0..i {
                    if l[(i + 1) * i / 2 + k] != 0. && l[(j + 1) * j / 2 + k] != 0. {
                        l[(j + 1) * j / 2 + i] -= l[(i + 1) * i / 2 + k] * l[(j + 1) * j / 2 + k]; 
                    }
                }
                l[(j + 1) * j / 2 + i] /= l[(i + 1) * i / 2 + i];
            }
        }
        let mut msg = Messenger::new("Solution of the system of equations", 1, 2 * n as i64, 1);
        // Вычисление y (Ly = b)
        let mut y: Array1<f64> = Array1::zeros(n);
        y[0] = self.b[0] / l[0];
        for i in 1..n {
            msg.add_progress();
            y[i] = self.b[i];
            for j in 0..i {
                y[i] -= y[j] * l[(i + 1) * i / 2 + j];
            }
            y[i] /= l[(i + 1) * i / 2 + i];
        }
        // Вычисление x (L(t)x = y)
        let mut x: Array1<f64> = Array1::zeros(n);
        x[n - 1] = y[n - 1] / l[l.len() - 1];
        for i in (0..n - 1).rev() {
            msg.add_progress();
            x[i] = y[i];
            for k in i + 1..n {
                x[i] -= x[k] * l[(k + 1) * k / 2 + i];
            }
            x[i] /= l[(i + 1) * i / 2 + i];
        }
        msg.stop();
        Ok(x)
    }
    // pub fn add_load(&mut self, index: usize, val: f64) -> Result<(), Error>  {
    //     if index >= self.nvtxs * self.blksze || index >= self.nvtxs * self.blksze  {
    //         return Err(Error::InvalidIndex);
    //     }
    //     self.b[index] += val;
    //     Ok(())
    // }
    // pub fn set_boundary_condition(&mut self, index: usize, val: f64) -> Result<(), Error> {
    //     if index >= self.nvtxs * self.blksze || index >= self.nvtxs * self.blksze  {
    //         return Err(Error::InvalidIndex);
    //     }
    //     for i in 0..self.nvtxs * self.blksze {
    //         if i != index {
    //             if self.a.set_value(index, i, 0.).is_err() {
    //                 continue;
    //             }
    //             if self.a.set_value(i, index, 0.).is_err() {
    //                 continue;
    //             }
    //         }
    //     }
    //     self.b[index] = val * self.a.get_value(index, index)?;
    //     Ok(())
    // }
}

impl Solver for LzhSolver {
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error> {
        self.a.add_value(index1, index2, value)
    }
    fn set_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error> {
        self.a.set_value(index1, index2, value)
    }
    fn set_vector_value(&mut self, index: usize, val: f64) -> Result<(), Error> {
        if index >= self.nvtxs * self.blksze {
            return Err(Error::InvalidIndex);
        }
        self.b[index] = val;
        Ok(())
    }
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), Error> {
        if index >= self.nvtxs * self.blksze {
            return Err(Error::InvalidIndex);
        }
        self.b[index] += val;
        Ok(())
    }
    fn set_result_value(&mut self, index: usize, val: f64) -> Result<(), Error> {
        self.a.clear_row(index)?;
        self.a.clear_col(index)?;
        self.a.set_value(index, index, 1.)?;
        self.b[index] = val;
        Ok(())
    }
}
