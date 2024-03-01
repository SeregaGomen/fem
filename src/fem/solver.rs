use super::error::FemError;
use super::sparse::{SparseMatrix, RussellSparseMatrix};
use super::mesh::Mesh;

pub trait FemSolver: Send {
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError>;
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError>; 
    fn mul_vector_value(&mut self, coef: f64);
    fn set_boundary_condition(&mut self, index: usize, val: f64) -> Result<(), FemError>;
    fn clear_matrix(&mut self);
    fn clear_vector(&mut self);
    fn solve(&mut self, eps: f64) -> Result<Vec<f64>, FemError>;
}

pub struct RussellSolver {
    matrix: RussellSparseMatrix,
    vector: Vec<f64>,
    boundary_condition: Vec<Option<f64>>,
}

impl RussellSolver {
    pub fn new(mesh: &Mesh) -> Result<Self, FemError> {
        let size = mesh.num_vertex * mesh.freedom;
        let nnz = 6 * mesh.nnz();
        Ok(Self { 
            matrix: RussellSparseMatrix::new(size, nnz)?, 
            vector: vec![0.0; size], 
            boundary_condition: vec![None; size],  
        })
    }
}

impl FemSolver for RussellSolver {
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError> {
        if self.boundary_condition[index1] == None && self.boundary_condition[index2] == None { self.matrix.add_value(index1, index2, value)? }
        Ok(())
    }
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.boundary_condition.len() { return Err(FemError::InvalidIndex) }
        if self.boundary_condition[index] == None { self.vector[index] += val }
        Ok(())
    }
    fn set_boundary_condition(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.boundary_condition.len() { return Err(FemError::InvalidIndex) }
        self.boundary_condition[index] = Some(val);
        Ok(())
    }
    fn mul_vector_value(&mut self, coef: f64) {
        for i in 0..self.vector.len() {
            self.vector[i] *= coef;
        }        
    }
    fn solve(&mut self, eps: f64) -> Result<Vec<f64>, FemError> {
        // Учет граничных условий
        for i in 0..self.boundary_condition.len() {
            if self.boundary_condition[i] != None {
                self.matrix.add_value(i, i, 1.0)?;
                self.vector[i] = self.boundary_condition[i].unwrap();
            }
        }
        self.matrix.solve(&self.vector, eps)
    }
    fn clear_matrix(&mut self) {
        self.matrix.clear();
    }
    fn clear_vector(&mut self) {
        self.vector.clear();
    }
}
