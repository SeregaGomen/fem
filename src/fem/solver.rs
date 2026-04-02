//use crate::fem::solver_russell::RussellSolver;
use crate::fem::sparse::SparseMatrix;
//use crate::fem::russell::RussellSparseMatrix;
use super::error::FemError;

pub trait FemSolver: Send {
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError>;
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError>; 
    fn mul_vector_value(&mut self, coef: f64);
    fn set_boundary_condition(&mut self, index: usize, val: f64) -> Result<(), FemError>;
    fn reset_matrix(&mut self);
    fn solve(&mut self, eps: f64) -> Result<Vec<f64>, FemError>;
}

pub struct GenericSolver<M: SparseMatrix> {
    matrix: M,
    vector: Vec<f64>,
    boundary_condition: Vec<Option<f64>>,
}

impl<M: SparseMatrix> GenericSolver<M> {
    pub fn new(matrix: M, size: usize) -> Self {
        Self {
            matrix,
            vector: vec![0.0; size],
            boundary_condition: vec![None; size],
        }
    }
}
impl<M: SparseMatrix + Send> FemSolver for GenericSolver<M> {
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError> {
        if self.boundary_condition[index1] == None && self.boundary_condition[index2] == None { self.matrix.add_value(index1, index2, value)? }
        Ok(())
    }
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.boundary_condition.len() { return Err(FemError::InvalidIndex) }
        if self.boundary_condition[index] == None { self.vector[index] += val }
        Ok(())
    }
    fn mul_vector_value(&mut self, coef: f64) {
        for i in 0..self.vector.len() {
            self.vector[i] *= coef;
        }
    }
    fn set_boundary_condition(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.boundary_condition.len() { return Err(FemError::InvalidIndex) }
        self.boundary_condition[index] = Some(val);
        Ok(())
    }
    fn reset_matrix(&mut self) {
        self.matrix.reset();
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
}

//pub type RussellSolver = GenericSolver<RussellSparseMatrix>;
