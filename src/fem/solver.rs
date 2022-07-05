use ndarray::Array1;
use super::error::FemError;
use super::sparse::{SparseMatrix, MapSparseMatrix, EnvSparseMatrix};
use super::mesh::Mesh;

pub trait FemSolver: Send {
    fn size(&self) -> usize;
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError>;
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError>; 
    fn solve(&mut self, eps: f64) -> Result<Array1<f64>, FemError>;
    fn mul_vector_value(&mut self, coef: f64);
    fn set_result_value(&mut self, index: usize, val: f64) -> Result<(), FemError>;
    fn clear_matrix(&mut self);
    fn clear_vector(&mut self);
}

pub struct LzhSolver {
    size: usize,
    a: MapSparseMatrix,
    b: Array1<f64>,
    bc: Vec<Option<f64>>,
}

pub struct EnvSolver {
    size: usize,
    a: EnvSparseMatrix,
    b: Array1<f64>,
    bc: Vec<Option<f64>>,
}

//unsafe impl Send for LzhSolver {}
//unsafe impl Send for EnvSolver {}

#[allow(dead_code)]
impl LzhSolver {
    pub fn new(mesh: &Mesh) -> Self {
        let size = mesh.num_vertex * mesh.freedom;
        Self { 
            size,
            a: MapSparseMatrix::new(mesh.num_vertex, mesh.freedom, &mesh.mesh_map), 
            b: Array1::zeros(size), 
            bc: vec![None; size], 
        }
    }
}

#[allow(dead_code)]
impl EnvSolver {
    pub fn new(mesh: &Mesh) -> Self {
        let size = mesh.num_vertex * mesh.freedom;
        Self { 
            size,
            a: EnvSparseMatrix::new(mesh.num_vertex, mesh.freedom, &mesh.mesh_map), 
            b: Array1::zeros(size), 
            bc: vec![None; size],  
        }
    }
}

impl FemSolver for LzhSolver {
    fn size(&self) -> usize {
        self.size
    }
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError> {
        if self.bc[index1] == None && self.bc[index2] == None { self.a.add_value(index1, index2, value)? }
        Ok(())
    }
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.size { return Err(FemError::InvalidIndex) }
        if self.bc[index] == None { self.b[index] += val }
        Ok(())
    }
    fn set_result_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.size() { return Err(FemError::InvalidIndex) }
        self.bc[index] = Some(val);
        Ok(())
    }
    fn mul_vector_value(&mut self, coef: f64) {
        for i in 0..self.b.len() {
            self.b[i] *= coef;
        }     
    }
    fn solve(&mut self, eps: f64) -> Result<Array1<f64>, FemError> {
        // Учет граничных условий
        for i in 0..self.size {
            if self.bc[i] != None {
                self.a.add_value(i, i, 1.0)?;
                self.b[i] = self.bc[i].unwrap();
            }
        }
        self.a.solve(&self.b, eps)
    }
    fn clear_matrix(&mut self) {
        self.a.clear();
    }
    fn clear_vector(&mut self) {
        self.b = Array1::zeros(self.b.len());
    }
}

impl FemSolver for EnvSolver {
    fn size(&self) -> usize {
        self.size
    }
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError> {
        if self.bc[index1] == None && self.bc[index2] == None { self.a.add_value(index1, index2, value)? }
        Ok(())
    }
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.size { return Err(FemError::InvalidIndex) }
        if self.bc[index] == None { self.b[index] += val }
        Ok(())
    }
    fn set_result_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.size() { return Err(FemError::InvalidIndex) }
        self.bc[index] = Some(val);
        Ok(())
    }
    fn mul_vector_value(&mut self, coef: f64) {
        for i in 0..self.b.len() {
            self.b[i] *= coef;
        }        
    }
    fn solve(&mut self, eps: f64) -> Result<Array1<f64>, FemError> {
        // Учет граничных условий
        for i in 0..self.size {
            if self.bc[i] != None {
                self.a.add_value(i, i, 1.0)?;
                self.b[i] = self.bc[i].unwrap();
            }
        }

        // println!();
        // for i in 0..self.size {
        //     for j in 0..self.size {
        //         print!("{:10.3} ", if self.a.get_value(i, j).is_err() { 0.0} else { self.a.get_value(i, j).unwrap() });
        //     }
        //     println!("{:10.3}", self.b[i]);
        // }
        // println!();


        self.a.solve(&self.b, eps)
    }
    fn clear_matrix(&mut self) {
        self.a.clear();
    }
    fn clear_vector(&mut self) {
        self.b = Array1::zeros(self.b.len());
    }
}
