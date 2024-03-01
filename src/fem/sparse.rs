use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry, LinSolKind};
use super::error::FemError;
use super::msg::Messenger;


pub trait SparseMatrix {
    fn add_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError>;
    fn solve(&mut self, rhs: &Vec<f64>, eps: f64) -> Result<Vec<f64>, FemError>;
    fn clear(&mut self);
}

// Russell 
pub struct RussellSparseMatrix {
    size: usize,
    data: SparseTriplet,
}

impl SparseMatrix for RussellSparseMatrix {
    fn add_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError> {
        self.data.put(index1, index2, value)?;
        Ok(())
    }
    fn solve(&mut self, rhs: &Vec<f64>, _: f64) -> Result<Vec<f64>, FemError> {
        // let mut solver = Solver::new(ConfigSolver::new())?;
        let mut solver = Solver::new(*ConfigSolver::new().lin_sol_kind(LinSolKind::Mmp))?;
        let mut msg = Messenger::new("Solution of the system of equations", 0, 0, 0);
        solver.initialize(&self.data)?;
        solver.factorize()?;
        let mut x = Vector::new(self.size);
        solver.solve(&mut x, & Vector::from(rhs))?;
        msg.stop();
        Ok(x.as_data().clone())
    }
    fn clear(&mut self) {
        self.data.reset();
    }
}

impl RussellSparseMatrix {
    pub fn new(size: usize, nnz: usize) -> Result<Self, FemError> {
        // Ok(Self{ size, data: SparseTriplet::new(size, size, nnz, Symmetry::General)? })
        Ok(Self{ size, data: SparseTriplet::new(size, size, nnz, Symmetry::No)? })
    }
}