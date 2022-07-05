use ndarray::Array1;
use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry, StrError};
use super::error::FemError;
use super::sparse::{SparseMatrix, MapSparseMatrix, EnvSparseMatrix};
use super::mesh::Mesh;

pub trait FemSolver: Send {
    fn size(&self) -> usize;
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError>;
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError>; 
    fn mul_vector_value(&mut self, coef: f64);
    fn set_result_value(&mut self, index: usize, val: f64) -> Result<(), FemError>;
    fn clear_matrix(&mut self);
    fn clear_vector(&mut self);
    fn solve(&mut self, eps: f64) -> Result<Vec<f64>, FemError>;
}

pub struct LzhSolver {
    size: usize,
    matrix: MapSparseMatrix,
    vector: Vec<f64>,
    boundary_condition: Vec<Option<f64>>,
}

pub struct EnvSolver {
    size: usize,
    matrix: EnvSparseMatrix,
    vector: Vec<f64>,
    boundary_condition: Vec<Option<f64>>,
}

pub struct RLSolver {
    size: usize,
    matrix: SparseTriplet,
    vector: Vec<f64>,
    boundary_condition: Vec<Option<f64>>,
}

#[allow(dead_code)]
impl LzhSolver {
    pub fn new(mesh: &Mesh) -> Self {
        let size = mesh.num_vertex * mesh.freedom;
        Self { 
            size,
            matrix: MapSparseMatrix::new(mesh.num_vertex, mesh.freedom, &mesh.mesh_map), 
            vector: vec![0.0; size], 
            boundary_condition: vec![None; size], 
        }
    }
}

#[allow(dead_code)]
impl EnvSolver {
    pub fn new(mesh: &Mesh) -> Self {
        let size = mesh.num_vertex * mesh.freedom;
        Self { 
            size,
            matrix: EnvSparseMatrix::new(mesh.num_vertex, mesh.freedom, &mesh.mesh_map), 
            vector: vec![0.0; size], 
            boundary_condition: vec![None; size],  
        }
    }
}

#[allow(dead_code)]
impl RLSolver {
    pub fn new(mesh: &Mesh) -> Self {
        let size = mesh.num_vertex * mesh.freedom;
        let nnz = 3 * mesh.nnz();
        Self { 
            size,
            matrix: SparseTriplet::new(size, size, nnz, Symmetry::General).unwrap(), 
            vector: vec![0.0; size], 
            boundary_condition: vec![None; size],  
        }
    }
}

impl FemSolver for LzhSolver {
    fn size(&self) -> usize {
        self.size
    }
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError> {
        if self.boundary_condition[index1] == None && self.boundary_condition[index2] == None { self.matrix.add_value(index1, index2, value)? }
        Ok(())
    }
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.size { return Err(FemError::InvalidIndex) }
        if self.boundary_condition[index] == None { self.vector[index] += val }
        Ok(())
    }
    fn set_result_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.size() { return Err(FemError::InvalidIndex) }
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
        for i in 0..self.size {
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
        self.vector = vec![0.0; self.size];
    }
}

impl FemSolver for EnvSolver {
    fn size(&self) -> usize {
        self.size
    }
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError> {
        if self.boundary_condition[index1] == None && self.boundary_condition[index2] == None { self.matrix.add_value(index1, index2, value)? }
        Ok(())
    }
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.size { return Err(FemError::InvalidIndex) }
        if self.boundary_condition[index] == None { self.vector[index] += val }
        Ok(())
    }
    fn set_result_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.size() { return Err(FemError::InvalidIndex) }
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
        for i in 0..self.size {
            if self.boundary_condition[i] != None {
                self.matrix.add_value(i, i, 1.0)?;
                self.vector[i] = self.boundary_condition[i].unwrap();
            }
        }

        // println!();
        // for i in 0..self.size {
        //     for j in 0..self.size {
        //         print!("{:10.3} ", if self.matrix.get_value(i, j).is_err() { 0.0} else { self.matrix.get_value(i, j).unwrap() });
        //     }
        //     println!("{:10.3}", self.vector[i]);
        // }
        // println!();


        self.matrix.solve(&self.vector, eps)
    }
    fn clear_matrix(&mut self) {
        self.matrix.clear();
    }
    fn clear_vector(&mut self) {
        self.vector = vec![0.0; self.size];
    }
}

impl FemSolver for RLSolver {
    fn size(&self) -> usize {
        self.size
    }
    fn add_matrix_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError> {
        if self.boundary_condition[index1] == None && self.boundary_condition[index2] == None { self.matrix.put(index1, index2, value)? }
        Ok(())
    }
    fn add_vector_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.size { return Err(FemError::InvalidIndex) }
        if self.boundary_condition[index] == None { self.vector[index] += val }
        Ok(())
    }
    fn set_result_value(&mut self, index: usize, val: f64) -> Result<(), FemError> {
        if index >= self.size() { return Err(FemError::InvalidIndex) }
        self.boundary_condition[index] = Some(val);
        Ok(())
    }
    fn mul_vector_value(&mut self, coef: f64) {
        for i in 0..self.size {
            self.vector[i] *= coef;
        }        
    }
    fn solve(&mut self, _: f64) -> Result<Vec<f64>, FemError> {
        // Учет граничных условий
        for i in 0..self.size {
            if self.boundary_condition[i] != None {
                self.matrix.put(i, i, 1.0)?;
                self.vector[i] = self.boundary_condition[i].unwrap();
            }
        }
        let mut solver = Solver::new(ConfigSolver::new()).unwrap();
        solver.initialize(&self.matrix).unwrap();
        solver.factorize().unwrap();
        let mut x = Vector::new(self.size);
        solver.solve(&mut x, & Vector::from(&self.vector)).unwrap();

        // let mut res = Array1::<f64>::zeros(self.size);
        // for i in 0..self.size {
        //     res[i] = x[i];
        // }
        let x = x.as_data().clone();
        Ok(x)
    }
    fn clear_matrix(&mut self) {
        self.matrix.reset();
    }
    fn clear_vector(&mut self) {
        self.vector = vec![0.0; self.size];
    }
}
