pub mod error;
mod fe;
mod util;
mod mesh;
mod sparse;
mod solver;
mod parser;
mod msg;

use rayon::prelude::*;
use std::sync::Mutex;
use bitflags::bitflags;
use ndarray::{Array1, Array2, prelude::*};
use std::time::Instant;
use error::FemError;
use mesh::Mesh;
// use solver::{Solver, LzhSolver, EnvSolver};
use solver::{Solver, LzhSolver};
// use solver::{Solver, EnvSolver};
use fe::FEType;
use parser::Parser;
use msg::Messenger;

bitflags! {
    pub struct Direct: u8 {
        const X = 0b00000001;
        const Y = 0b00000010;
        const Z = 0b00000100;
    }
}

#[derive(PartialEq)]
enum ParamType {
    BoundaryCondition,  // Граничное условие, заданноe выражением
    VolumeLoad,         // Объемная нагрузка
    SurfaceLoad,        // Поверхностная ...
    ConcentratedLoad,   // Сосредоточенная ...
    PressureLoad,       // Нагрузка давлением
}

struct Parameter<'a> {
    p_type: ParamType,
    value: &'a str,
    predicate: &'a str,
    value_fun: Option<fn(f64, f64, f64) -> f64>,
    predicate_fun: Option<fn(f64, f64, f64) -> bool>,
    direct: Direct,
}

impl<'a> Parameter<'a> {
    fn new_str(p_type: ParamType, value: &'a str, predicate: &'a str, direct: Direct) -> Self {
        Self { p_type, value, predicate, value_fun: None, predicate_fun: None, direct }
    }
    fn new_fun(p_type: ParamType, value_f: fn(f64, f64, f64) -> f64, predicate_f: fn(f64, f64, f64) -> bool, direct: Direct) -> Self {
        Self { p_type, value: "", predicate: "", value_fun: Some(value_f), predicate_fun: Some(predicate_f), direct }
    }
    fn get_value(&self, x: &ArrayView1<f64>) -> Result<f64, FemError> {
        match self.value_fun {
            None => {
                let var_name = array!["x", "y", "shell-tube-3z"];
                let mut parser = Parser::new();
                for i in 0..x.len() {
                    parser.set_variable(var_name[i], x[i]);    
                }
                parser.set_expression(&self.value)?;
                parser.value()
            }
            Some(fun) => Ok(fun(x[0], if x.len() > 1 { x[1] } else { 0. }, if x.len() > 2 { x[2] } else { 0. } )),
        }
    }
    fn get_predicate(&self, x: &ArrayView1<f64>) -> Result<bool, FemError> {
        match self.predicate_fun {
            None => {
                if self.predicate.len() == 0 {
                    Ok(true)
                } else {
                    let var_name = array!["x", "y", "z"];
                    let mut parser = Parser::new();
                    for i in 0..x.len() {
                        parser.set_variable(var_name[i], x[i]);    
                    }
                    parser.set_expression(&self.predicate)?;
                    if parser.value()? == 1. { Ok(true) } else { Ok(false) }
                }
            }
            Some(fun) => Ok(fun(x[0], if x.len() > 1 { x[1] } else { 0. }, if x.len() > 2 { x[2] } else { 0. } )),
        }
    }
}

struct FEMParameter<'a> {
    param: Vec<Parameter<'a>>,
    e: [f64; 2],
    thk: f64,
    eps: f64,
    nthreads: usize,
}

impl<'a> FEMParameter<'a> {
    fn new() -> Self {
        Self { param: Vec::new(), e: [0., 0.], thk: 0., eps: 1.0e-6, nthreads: 1 }
    }
    fn find_parameter(&self, param: ParamType) -> bool {
        for i in &self.param {
            if i.p_type == param {
                return true;
            }
        }
        false
    }
}

pub struct FEM<'a> {
    mesh: Mesh,
    param: FEMParameter<'a>,
}

#[allow(dead_code)]
impl<'a> FEM<'a> {
    pub fn new(mesh_name: &str) -> Result<Self, FemError> {
        let mesh = Mesh::new(mesh_name)?;
        Ok(Self {mesh, param: FEMParameter::new()})
    }
    pub fn set_num_threads(&mut self, nthreads: usize) {
        self.param.nthreads = nthreads;
    }
    pub fn set_young_modulus(&mut self, e: f64) {
        self.param.e[0] = e;
    }
    pub fn set_poisons_ratio(&mut self, m: f64) {
        self.param.e[1] = m;
    }
    pub fn set_thickness(&mut self, thk: f64) {
        self.param.thk = thk;
    }
    pub fn set_eps(&mut self, eps: f64) {
        self.param.eps = eps;
    }
    pub fn add_boundary_condition(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.param.push(Parameter::new_str(ParamType::BoundaryCondition, value, predicate, direct));
    }
    pub fn add_boundary_condition_fun(&mut self, value_fun: fn(f64, f64, f64) -> f64, predicate_fun: fn(f64, f64, f64) -> bool, direct: Direct) {
        self.param.param.push(Parameter::new_fun(ParamType::BoundaryCondition, value_fun, predicate_fun, direct));
    }
    pub fn add_concentrated_load(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.param.push(Parameter::new_str(ParamType::ConcentratedLoad, value, predicate, direct));
    }
    pub fn add_volume_load(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.param.push(Parameter::new_str(ParamType::VolumeLoad, value, predicate, direct));
    }
    pub fn add_volume_load_fun(&mut self, value_fun: fn(f64, f64, f64) -> f64, predicate_fun: fn(f64, f64, f64) -> bool, direct: Direct) {
        self.param.param.push(Parameter::new_fun(ParamType::VolumeLoad, value_fun, predicate_fun, direct));
    }
    pub fn add_surface_load(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.param.push(Parameter::new_str(ParamType::SurfaceLoad, value, predicate, direct));
    }
    pub fn add_surface_load_fun(&mut self, value_fun: fn(f64, f64, f64) -> f64, predicate_fun: fn(f64, f64, f64) -> bool, direct: Direct) {
        self.param.param.push(Parameter::new_fun(ParamType::SurfaceLoad, value_fun, predicate_fun, direct));
    }
    pub fn add_pressure_load(&mut self, value: &'a str, predicate: &'a str) {
        self.param.param.push(Parameter::new_str(ParamType::PressureLoad, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_pressure_load_fun(&mut self, value_fun: fn(f64, f64, f64) -> f64, predicate_fun: fn(f64, f64, f64) -> bool) {
        self.param.param.push(Parameter::new_fun(ParamType::PressureLoad, value_fun, predicate_fun, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn is_2d(&self) -> bool {
        self.mesh.is_2d()
    }
    pub fn is_3d(&self) -> bool {
        self.mesh.is_3d()
    }
    pub fn generate(&mut self, res_name: &str) -> Result<(), FemError> {
        rayon::ThreadPoolBuilder::new().num_threads(self.param.nthreads).build_global().unwrap();
        let time = Instant::now();
        let mut solver = Mutex::new(LzhSolver::new(&self.mesh));
        // let mut solver = Mutex::new(EnvSolver::new(&self.mesh));
        self.set_global_matrix(&mut solver)?;
        self.set_concentrated_load(&mut solver)?;
        self.set_volume_load(&mut solver)?;
        self.set_surface_load(&mut solver)?;
        self.set_boundary_condition(&mut solver)?;
        self.calc_results(&solver.lock().unwrap().solve(self.param.eps)?, res_name)?;
        println!("Lead time: {:.2?}", time.elapsed());
        Ok(())
    }
    fn set_global_matrix(&self, solver: &mut Mutex<impl Solver>) -> Result<(), FemError> {
        let msg = Mutex::new(Messenger::new("Generate global stiffness matrix", 1, self.mesh.num_fe as i64, 5));
        (0..self.mesh.num_fe).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
            msg.lock().unwrap().add_progress();
            let fe = self.calc_fe_matrix(i)?;
            for j in 0..fe.shape()[0] {
                for k in j..fe.shape()[1] {
                    let index1 = self.mesh.fe[[i, j / self.mesh.freedom]] * self.mesh.freedom + j % self.mesh.freedom;
                    let index2 = self.mesh.fe[[i, k / self.mesh.freedom]] * self.mesh.freedom + k % self.mesh.freedom;
                    solver.lock().unwrap().add_matrix_value(index1, index2, fe[[j, k]])?;
                    if j != k {
                        solver.lock().unwrap().add_matrix_value(index2, index1, fe[[j, k]])?;
                    }
                }
            }
            Ok(())
        })?;
        msg.lock().unwrap().stop();
        Ok(())
    }
    fn num_results(&self) -> usize {
        match self.mesh.fe_type {
            FEType::FE1D2 => 3,
            FEType::FE2D3 | FEType::FE2D4 => 8,
            FEType::FE3D4 | FEType::FE3D8 => 15,
            FEType::FE2D3S | FEType::FE2D4S => 18,
        }
    }
    fn fun_names(&self) -> Vec<&str> {
        match self.mesh.fe_type {
            FEType::FE1D2 => vec![ "U", "Exx", "Sxx" ],
            FEType::FE2D3 | FEType::FE2D4 => vec![ "U", "V", "Exx", "Eyy", "Exy", "Sxx", "Syy", "Sxy" ],
            FEType::FE3D4 | FEType::FE3D8 => vec![ "U", "V", "W", "Exx", "Eyy", "Ezz", "Exy", "Exz", "Eyz", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz" ],
            FEType::FE2D3S | FEType::FE2D4S => vec![ "U", "V", "W", "Tx", "Ty", "Tz", "Exx", "Eyy", "Ezz", "Exy", "Exz", "Eyz", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz" ],
        }
    }
    fn print_summary(&self, res: &Array2<f64>) {
        let names = self.fun_names();
        println!("Fun:\tmin\t\tmax");
        for i in 0..res.shape()[0] {
            // println!("{}:\t{:+.6e}\t{:+.6e}", names[i], min_max[[i, 0]], min_max[[i, 1]]);
            println!("{}:\t{}\t{}", names[i], util::fmt_f64(util::get_min(res.row(i)), 12, 5, 2), util::fmt_f64(util::get_max(res.row(i)), 12, 5, 2));
        }
    }
    fn calc_results(&self, u: &Array1<f64>, res_name: &str) -> Result<(), FemError> {
        let fe_size = self.mesh.fe.shape()[1];
        let freedom = self.mesh.freedom;
        let res = Mutex::new(Array2::zeros((self.num_results(), self.mesh.num_vertex)));
        let counter = Mutex::new(Array1::<i32>::zeros(self.mesh.num_vertex));
        // Копирование результатов расчета (перемещений)
        for i in 0..freedom {
            for j in 0..self.mesh.num_vertex {
                res.lock().unwrap()[[i, j]] = u[j * freedom + i];
            }
        }
        // Вычисление деформаций, напряжений, ...
        let msg = Mutex::new(Messenger::new("Calculation of standard finite element results", 1, self.mesh.num_fe as i64, 5));
        (0..self.mesh.num_fe).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
            msg.lock().unwrap().add_progress();
            // Формируем вектор перемещений для текущего КЭ
            let mut fe_u = Array1::zeros(fe_size * freedom);
            for j in 0..fe_size {
                for k in 0..freedom {
                    fe_u[j * freedom + k] = u[freedom * self.mesh.fe[[i, j]] + k];
                }
            }
            let fe_res = self.calc_fe_results(i, fe_u)?; 
            for j in 0..self.num_results() - freedom {
                for k in 0..fe_size {
                    res.lock().unwrap()[[j + freedom, self.mesh.fe[[i, k]]]] += fe_res[[j, k]];
                    if j == 0 {
                        counter.lock().unwrap()[self.mesh.fe[[i, k]]] += 1;
                    }
                }
            }
            Ok(())
        })?;
        let mut res = res.lock().unwrap();
        let counter = counter.lock().unwrap();
        // Осредняем результаты
        for i in freedom..self.num_results() {
            for j in 0..self.mesh.num_vertex {
                res[[i, j]] /= counter[j] as f64;
            }
        }
        self.print_summary(&res);
        self.save_results(&res, res_name)
    }
    fn set_boundary_condition(&mut self, solver: &mut Mutex<impl Solver>) -> Result<(), FemError> {
        let msg = Mutex::new(Messenger::new("Using of boundary conditions", 1, self.mesh.num_vertex as i64, 5));
        (0..self.mesh.num_vertex).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
            msg.lock().unwrap().add_progress();
            for it in &self.param.param {
                if it.p_type == ParamType::BoundaryCondition {
                    let x = self.mesh.x.row(i);
                    if it.get_predicate(&x)? == false {
                        continue;
                    }
                    let val = it.get_value(&x)?;
                    if it.direct.contains(Direct::X) {
                        // solver.lock().unwrap().set_result_value((i + l) * self.mesh.freedom + 0, val)?;    
                        solver.lock().unwrap().set_result_value((i) * self.mesh.freedom + 0, val)?;    
                    }
                    if it.direct.contains(Direct::Y) && (self.mesh.is_2d() || self.mesh.is_3d()) {
                        solver.lock().unwrap().set_result_value((i) * self.mesh.freedom + 1, val)?;    
                    }
                    if it.direct.contains(Direct::Z) && self.mesh.is_3d() {
                        solver.lock().unwrap().set_result_value((i) * self.mesh.freedom + 2, val)?;    
                    }
                }
            }
            Ok(())    
        })?;
        msg.lock().unwrap().stop();
        Ok(())
    }
    fn set_concentrated_load(&mut self, solver: &mut Mutex<impl Solver>) -> Result<(), FemError> {
        if self.param.find_parameter(ParamType::ConcentratedLoad) {
            let msg = Mutex::new(Messenger::new("Calculation of concentrated loads", 1, self.mesh.num_vertex as i64, 5));
            (0..self.mesh.num_vertex).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
                msg.lock().unwrap().add_progress();
                for it in &self.param.param {
                    if it.p_type == ParamType::ConcentratedLoad {
                        let x = self.mesh.x.row(i);
                        if it.get_predicate(&x)? == false {
                            continue;
                        }
                        let val = it.get_value(&x)?;
                        if it.direct.contains(Direct::X) {
                            solver.lock().unwrap().set_vector_value((i) * self.mesh.freedom + 0, val)?;    
                        }
                        if it.direct.contains(Direct::Y) && (self.mesh.is_2d() || self.mesh.is_3d()) {
                            solver.lock().unwrap().set_vector_value((i) * self.mesh.freedom + 1, val)?;    
                        }
                        if it.direct.contains(Direct::Z) && self.mesh.is_3d() {
                            solver.lock().unwrap().set_vector_value((i) * self.mesh.freedom + 2, val)?;    
                        }
                    }
                }
                Ok(())
            })?;
        }
        Ok(())
    }
    fn set_volume_load(&mut self, solver: &mut Mutex<impl Solver>) -> Result<(), FemError> {
        if self.param.find_parameter(ParamType::VolumeLoad) {
            let msg = Mutex::new(Messenger::new("Calculation of volume loads", 1, self.mesh.num_fe as i64, 5));
            (0..self.mesh.num_fe).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
                msg.lock().unwrap().add_progress();
                for it in &self.param.param {
                    if it.p_type == ParamType::VolumeLoad {
                        let x = self.mesh.get_fe_coord(i);
                        if !self.check_elem(&x, &it)? {
                            continue;
                        }
                        let share = self.volume_load_share() * self.mesh.fe_volume(i); 
                        for j in 0..self.mesh.fe.shape()[1] {
                            let val = it.get_value(&x.row(j))? * share[j];
                            if it.direct.contains(Direct::X) {
                                solver.lock().unwrap().add_vector_value(self.mesh.fe[[i, j]] * self.mesh.freedom + 0, val)?;    
                            }
                            if it.direct.contains(Direct::Y) && (self.mesh.is_2d() || self.mesh.is_3d()) {
                                solver.lock().unwrap().add_vector_value(self.mesh.fe[[i, j]] * self.mesh.freedom + 1, val)?;    
                            }
                            if it.direct.contains(Direct::Z) && self.mesh.is_3d() {
                                solver.lock().unwrap().add_vector_value(self.mesh.fe[[i, j]] * self.mesh.freedom + 2, val)?;    
                            }
                        }
                    }
                }
                Ok(())
            })?;
        }
        Ok(())
    }
    fn set_surface_load(&mut self, solver: &mut Mutex<impl Solver>) -> Result<(), FemError> {
        if self.param.find_parameter(ParamType::PressureLoad) || self.param.find_parameter(ParamType::SurfaceLoad) {
            let msg = Mutex::new(Messenger::new("Calculation of surface loads", 1, self.mesh.num_be as i64, 5));
            (0..self.mesh.num_be).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
                msg.lock().unwrap().add_progress();
                for it in &self.param.param {
                    if it.p_type == ParamType::PressureLoad || it.p_type == ParamType::SurfaceLoad {
                        let x = self.mesh.get_be_coord(i);
                        if !self.check_elem(&x, &it)? {
                            continue;
                        }
                        let share = self.surface_load_share() * self.mesh.be_volume(i); 
                        let normal = if it.p_type == ParamType::PressureLoad { self.mesh.be_normal(i) } else { array![1., 1., 1.] };
                        for j in 0..self.mesh.be.shape()[1] {
                            let val = it.get_value(&x.row(j))? * share[j];
                            if it.direct.contains(Direct::X) {
                                solver.lock().unwrap().add_vector_value(self.mesh.be[[i, j]] * self.mesh.freedom + 0, val * normal[0])?;    
                            }
                            if it.direct.contains(Direct::Y) && (self.mesh.is_2d() || self.mesh.is_3d()) {
                                solver.lock().unwrap().add_vector_value(self.mesh.be[[i, j]] * self.mesh.freedom + 1, val * normal[1])?;    
                            }
                            if it.direct.contains(Direct::Z) && self.mesh.is_3d() {
                                solver.lock().unwrap().add_vector_value(self.mesh.be[[i, j]] * self.mesh.freedom + 2, val * normal[2])?;    
                            }
                        }
                    }
                }
                Ok(())
            })?;
        }
        Ok(())
    }
    fn check_elem(&self, x: &Array2<f64>, param: &Parameter) -> Result<bool, FemError> {
        for i in 0..x.shape()[0] {
            if param.get_predicate(&x.row(i))? == false {
                return Ok(false);
            }
        }
        Ok(true)
    }
    fn volume_load_share(&self) -> Array1<f64> {
        match self.mesh.fe_type {   
            FEType::FE1D2 => array![ 0.5, 0.5 ],
            FEType::FE2D3 | FEType::FE2D3S => array![ 0.33333333333, 0.33333333333, 0.33333333333 ],
            FEType::FE2D4 | FEType::FE2D4S | FEType::FE3D4 => array![ 0.25, 0.25, 0.25, 0.25 ],
            FEType::FE3D8 => array![ 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125 ],
        } 
    }
    fn surface_load_share(&self) -> Array1<f64> {
        match self.mesh.fe_type {   
            FEType::FE1D2 => array![ 1.0 ],
            FEType::FE2D3 | FEType::FE2D4 => array![ 0.5, 0.5 ],
            FEType::FE2D3S | FEType::FE3D4 => array![ 0.333333333333, 0.333333333333, 0.333333333333 ],
            FEType::FE2D4S | FEType::FE3D8 => array![ 0.25, 0.25, 0.25, 0.25 ],
        } 
    }
    fn calc_fe_matrix(&self, i: usize) -> Result<Array2<f64>, FemError> {
        use crate::fem::fe::fe1d::FiniteElement1D;
        use crate::fem::fe::fe2d::FiniteElement2D;
        use crate::fem::fe::fe3d::FiniteElement3D;

        match self.mesh.fe_type {   
            FEType::FE1D2 => {
                let mut fe = fe::fe1d::FE1D2::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i));
                fe.generate()
            }
            FEType::FE2D3 | FEType::FE2D3S  => {
                let mut fe = fe::fe2d::FE2D3::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i), if self.mesh.fe_type == FEType::FE2D3 { false } else { true });
                fe.generate()
            }
            FEType::FE2D4 | FEType::FE2D4S => {
                let mut fe = fe::fe2d::FE2D4::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i), if self.mesh.fe_type == FEType::FE2D4 { false } else { true });
                fe.generate()
            }
            FEType::FE3D4 => {
                let mut fe = fe::fe3d::FE3D4::new(self.param.e, self.mesh.get_fe_coord(i));
                fe.generate()
            }
            FEType::FE3D8 => {
                let mut fe = fe::fe3d::FE3D8::new(self.param.e, self.mesh.get_fe_coord(i));
                fe.generate()
            }
        } 
    }
    fn calc_fe_results(&self, i: usize, u: Array1<f64>) -> Result<Array2<f64>, FemError> {
        use crate::fem::fe::fe1d::FiniteElement1D;
        use crate::fem::fe::fe2d::FiniteElement2D;
        use crate::fem::fe::fe3d::FiniteElement3D;

        match self.mesh.fe_type {   
            FEType::FE1D2 => {
                let fe = fe::fe1d::FE1D2::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i));
                fe.calc(&u)
            }
            FEType::FE2D3 | FEType::FE2D3S => {
                let fe = fe::fe2d::FE2D3::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i), if self.mesh.fe_type == FEType::FE2D3 { false } else { true });
                fe.calc(&u)
            }
            FEType::FE2D4 | FEType::FE2D4S => {
                let fe = fe::fe2d::FE2D4::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i), if self.mesh.fe_type == FEType::FE2D4 { false } else { true });
                fe.calc(&u)
            }
            FEType::FE3D4 => {
                let fe = fe::fe3d::FE3D4::new(self.param.e, self.mesh.get_fe_coord(i));
                fe.calc(&u)
            }
            FEType::FE3D8 => {
                let fe = fe::fe3d::FE3D8::new(self.param.e, self.mesh.get_fe_coord(i));
                fe.calc(&u)
            }
        } 
    }
    fn save_results(&self, res: &Array2<f64>, file_name: &str) -> Result<(), FemError> {
        use chrono::{Datelike, Timelike, Utc};
        use std::fs::File;
        use std::io::{BufWriter, prelude::*};
                
        let file: File = match File::create(file_name) {
            Err(_) => return Err(FemError::OpenFile),
            Ok(file) => file,
        };
        let mut stream = BufWriter::new(file);
        // Запись сигнатуры
        if !write!(stream, "FEM Solver Results File\n").is_ok() {
            return Err(FemError::WriteFile);
        }
        // Вывод сетки
        if !write!(stream, "Mesh\n").is_ok() {
            return Err(FemError::WriteFile);
        }
        if !write!(stream, "{}\n", self.mesh.fe_type).is_ok() {
            return Err(FemError::WriteFile);
        }
        if !write!(stream, "{}\n", self.mesh.num_vertex).is_ok() {
            return Err(FemError::WriteFile);
        }
        for i in 0..self.mesh.x.shape()[0] {
            for j in 0..self.mesh.x.shape()[1] {
                if !write!(stream, "{} ", self.mesh.x[[i, j]]).is_ok() {
                    return Err(FemError::WriteFile);
                }
            }
            if !write!(stream, "\n").is_ok() {
                return Err(FemError::WriteFile);
            }
        }
        if !write!(stream, "{}\n", self.mesh.num_fe).is_ok() {
            return Err(FemError::WriteFile);
        }
        for i in 0..self.mesh.fe.shape()[0] {
            for j in 0..self.mesh.fe.shape()[1] {
                if !write!(stream, "{} ", self.mesh.fe[[i, j]]).is_ok() {
                    return Err(FemError::WriteFile);
                }
            }
            if !write!(stream, "\n").is_ok() {
                return Err(FemError::WriteFile);
            }
        }
        if self.mesh.is_shell() {
            if !write!(stream, "0\n").is_ok() {
                return Err(FemError::WriteFile);
            }
        } else {
            if !write!(stream, "{}\n", self.mesh.num_be).is_ok() {
                return Err(FemError::WriteFile);
            }
            for i in 0..self.mesh.be.shape()[0] {
                for j in 0..self.mesh.be.shape()[1] {
                    if !write!(stream, "{} ", self.mesh.be[[i, j]]).is_ok() {
                        return Err(FemError::WriteFile);
                    }
                }
                if !write!(stream, "\n").is_ok() {
                    return Err(FemError::WriteFile);
                }
            }
        }
        // Запись результатов расчета
        if !write!(stream, "Results\n").is_ok() {
            return Err(FemError::WriteFile);
        }
        let now = Utc::now();
        if !write!(stream, "{:02}.{:02}.{:4} - {:02}:{:02}:{:02}\n",now.day(), now.month(), now.year(), now.hour(), now.minute(), now.second()).is_ok() {
            return Err(FemError::WriteFile);
        }
        if !write!(stream, "{}\n", self.num_results()).is_ok() {
            return Err(FemError::WriteFile);
        }
        for i in 0..res.shape()[0] {
            if !write!(stream, "{}\n", self.fun_names()[i]).is_ok() {
                return Err(FemError::WriteFile);
            }
            if !write!(stream, "0\n{}\n", res.shape()[1]).is_ok() {
                return Err(FemError::WriteFile);
            }
            for j in 0..res.shape()[1] {
                if !write!(stream, "{}\n", res[[i, j]]).is_ok() {
                    return Err(FemError::WriteFile);
                }
            }
        }
        Ok(())
    }
}