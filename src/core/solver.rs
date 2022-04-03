use ndarray::{Array1};
use crate::error::{ErrorCode, Error, error};
use super::sparse::{SparseMatrix, MapSparseMatrix, EnvSparseMatrix};
use super::mesh::Mesh;

pub trait Solver: Send {
    fn size(&self) -> usize;
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error>;
    fn set_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error>;
    fn get_matrix_value(&self, index1: usize, index2: usize) -> Result<f64, Error>;
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), Error>; 
    fn set_vector_value(&mut self, index: usize, val: f64) -> Result<(), Error>;
    fn solve(&mut self, eps: f64) -> Result<Array1<f64>, Error>;
    fn set_result_value(&mut self, index: usize, val: f64) -> Result<(), Error> {
        for i in 0..self.size() {
            if i != index {
                if self.set_matrix_value(index, i, 0.).is_err() {
                    continue;
                }
                self.set_matrix_value(i, index, 0.)?;
            }
        }
        self.set_vector_value(index, val * self.get_matrix_value(index, index)?)?;
        Ok(())
    }
}

pub struct LzhSolver {
    size: usize,
    a: MapSparseMatrix,
    b: Array1<f64>,
}

pub struct EnvSolver {
    size: usize,
    a: EnvSparseMatrix,
    b: Array1<f64>,
}

//unsafe impl Send for LzhSolver {}
//unsafe impl Send for EnvSolver {}

#[allow(dead_code)]
impl LzhSolver {
    pub fn new(mesh: &Mesh) -> Self {
        Self { 
            size: mesh.num_vertex * mesh.freedom,
            a: MapSparseMatrix::new(mesh.num_vertex, mesh.freedom, &mesh.mesh_map), 
            b: Array1::zeros(mesh.num_vertex * mesh.freedom), 
        }
    }
}

#[allow(dead_code)]
impl EnvSolver {
    pub fn new(mesh: &Mesh) -> Self {
        Self { 
            size: mesh.num_vertex * mesh.freedom,
            a: EnvSparseMatrix::new(mesh.num_vertex, mesh.freedom, &mesh.mesh_map), 
            b: Array1::zeros(mesh.num_vertex * mesh.freedom), 
        }
    }
}

impl Solver for LzhSolver {
    fn size(&self) -> usize {
        self.size
    }
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error> {
        self.a.add_value(index1, index2, value)
    }
    fn set_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error> {
        self.a.set_value(index1, index2, value)
    }
    fn get_matrix_value(&self, index1: usize, index2: usize) -> Result<f64, Error> {
        self.a.get_value(index1, index2)
    }
    fn set_vector_value(&mut self, index: usize, val: f64) -> Result<(), Error> {
        if index >= self.size {
            return Err(error(ErrorCode::InvalidIndex));
        }
        self.b[index] = val;
        Ok(())
    }
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), Error> {
        if index >= self.size {
            return Err(error(ErrorCode::InvalidIndex));
        }
        self.b[index] += val;
        Ok(())
    }
    fn solve(&mut self, eps: f64) -> Result<Array1<f64>, Error> {
        self.a.solve(&self.b, eps)
    }
}

impl Solver for EnvSolver {
    fn size(&self) -> usize {
        self.size
    }
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error> {
        self.a.add_value(index1, index2, value)
    }
    fn set_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), Error> {
        self.a.set_value(index1, index2, value)
    }
    fn get_matrix_value(&self, index1: usize, index2: usize) -> Result<f64, Error> {
        self.a.get_value(index1, index2)
    }
    fn set_vector_value(&mut self, index: usize, val: f64) -> Result<(), Error> {
        if index >= self.size {
            return Err(error(ErrorCode::InvalidIndex));
        }
        self.b[index] = val;
        Ok(())
    }
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), Error> {
        if index >= self.size {
            return Err(error(ErrorCode::InvalidIndex));
        }
        self.b[index] += val;
        Ok(())
    }
    fn solve(&mut self, eps: f64) -> Result<Array1<f64>, Error> {
        self.a.solve(&self.b, eps)
    }
}
