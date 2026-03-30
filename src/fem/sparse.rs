use super::error::FemError;


pub trait SparseMatrix {
    fn add_value(&mut self, index1: usize, index2: usize, value: f64) -> Result<(), FemError>;
    fn solve(&mut self, rhs: &Vec<f64>, eps: f64) -> Result<Vec<f64>, FemError>;
    fn clear(&mut self);
    fn reset(&mut self);
}

