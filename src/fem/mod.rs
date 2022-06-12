pub mod error;
pub mod mesh;
mod fe;
mod util;
mod sparse;
mod solver;
mod parser;
mod msg;

use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use rayon::prelude::*;
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

enum ParamValueType {
    String,
    Function,
}

struct ParamValue<'a> {
    value_type: ParamValueType,
    value_str: &'a str,
    value_fun: Option<fn(f64, f64, f64) -> f64>,
}

impl<'a> ParamValue<'a> {
    fn new_str(value: &'a str) -> Self {
        Self { value_type: ParamValueType::String, value_str: value, value_fun: None }
    }
    fn new_fun(value: fn(f64, f64, f64) -> f64) -> Self {
        Self { value_type: ParamValueType::Function, value_str: "", value_fun: Some(value) }
    }
} 

#[derive(PartialEq)]
pub enum ParamType {
    BoundaryCondition,  // Граничное условие, заданноe выражением
    VolumeLoad,         // Объемная нагрузка
    SurfaceLoad,        // Поверхностная ...
    ConcentratedLoad,   // Сосредоточенная ...
    PressureLoad,       // Нагрузка давлением
    Thickness,          // Толщина (сечение) элемента
    YoungModulus,       // Модуль Юнга
    PoissonsRatio,      // Коэффициент Пуассона
    StressStrainCurve,  // Диаграмма материала
}

pub struct Parameter<'a> {
    p_type: ParamType,
    value: ParamValue<'a>,
    predicate: ParamValue<'a>,
    direct: Direct,
}

impl<'a> Parameter<'a> {
    fn new_str(p_type: ParamType, val: &'a str, pred: &'a str, direct: Direct) -> Self {
        let value = ParamValue::new_str(val);
        let predicate = ParamValue::new_str(pred);
        Self { p_type, value, predicate, direct }
    }
    fn new_fun(p_type: ParamType, val: fn(f64, f64, f64) -> f64, pred: fn(f64, f64, f64) -> f64, direct: Direct) -> Self {
        let value = ParamValue::new_fun(val);
        let predicate = ParamValue::new_fun(pred);
        Self { p_type, value, predicate, direct }
    }
    fn get_value(&self, x: &ArrayView1<f64>, variables: &HashMap<&'a str, f64>) -> Result<f64, FemError> {
        match self.value.value_type {
            ParamValueType::String => {
                let var_name = ["x", "y", "z"];
                let mut parser = Parser::new();
                for i in 0..x.len() {
                    parser.set_variable(var_name[i], x[i]);    
                }
                for i in variables {
                    parser.set_variable(i.0, *i.1);    
                }
                parser.set_expression(self.value.value_str)?;
                parser.value()
            }
            ParamValueType::Function => {
                if let Some(fun) = self.value.value_fun  {
                    Ok(fun(x[0], if x.len() > 1 { x[1] } else { 0. }, if x.len() > 2 { x[2] } else { 0. })) 
                } else {
                    Err(FemError::InternalError)
                }
            }
        }
    }
    fn get_predicate(&self, x: &ArrayView1<f64>, variables: &HashMap<&'a str, f64>) -> Result<bool, FemError> {
        match self.predicate.value_type {
            ParamValueType::String => {
                if self.predicate.value_str.len() == 0 {
                    Ok(true)
                } else {
                    let var_name = ["x", "y", "z"];
                    let mut parser = Parser::new();
                    for i in 0..x.len() {
                        parser.set_variable(var_name[i], x[i]);    
                    }
                    for i in variables {
                        parser.set_variable(i.0, *i.1);    
                    }
                    parser.set_expression(&self.predicate.value_str)?;
                    Ok(if parser.value()? == 1. { true } else { false })
                }
            }
            ParamValueType::Function => {
                let res = if let Some(fun) = self.value.value_fun {
                    fun(x[0], if x.len() > 1 { x[1] } else { 0. }, if x.len() > 2 { x[2] } else { 0. } )
                } else {
                    return Err(FemError::InternalError)
                };
                Ok(if res == 1. {true} else {false})
            }
        }
    }
}

pub struct FEMParameter<'a> {
    param: Vec<Parameter<'a>>,
    eps: f64,
    nthreads: usize,
    variables: HashMap<&'a str, f64>,
}

impl<'a> FEMParameter<'a> {
    pub fn new() -> Self {
        Self { param: Vec::new(), eps: 1.0e-6, nthreads: 1, variables: HashMap::new() }
    }
    fn find_parameter(&self, param: ParamType) -> bool {
        for i in &self.param {
            if i.p_type == param {
                return true;
            }
        }
        false
    }
    pub fn set_num_threads(&mut self, nthreads: usize) {
        self.nthreads = nthreads;
    }
    pub fn add_young_modulus_str(&mut self, value: &'a str, predicate: &'a str) {
        self.param.push(Parameter::new_str(ParamType::YoungModulus, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_young_modulus_fun(&mut self, value: fn(f64, f64, f64) -> f64, predicate: fn(f64, f64, f64) -> f64) {
        self.param.push(Parameter::new_fun(ParamType::YoungModulus, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_poisons_ratio_str(&mut self, value: &'a str, predicate: &'a str) {
        self.param.push(Parameter::new_str(ParamType::PoissonsRatio, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_poisons_ratio_fun(&mut self, value: fn(f64, f64, f64) -> f64, predicate: fn(f64, f64, f64) -> f64) {
        self.param.push(Parameter::new_fun(ParamType::PoissonsRatio, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_thickness_str(&mut self, value: &'a str, predicate: &'a str) {
        self.param.push(Parameter::new_str(ParamType::Thickness, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_thickness_fun(&mut self, value: fn(f64, f64, f64) -> f64, predicate: fn(f64, f64, f64) -> f64) {
        self.param.push(Parameter::new_fun(ParamType::Thickness, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn set_eps(&mut self, eps: f64) {
        self.eps = eps;
    }
    pub fn add_boundary_condition_str(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.push(Parameter::new_str(ParamType::BoundaryCondition, value, predicate, direct));
    }
    pub fn add_boundary_condition_fun(&mut self, value: fn(f64, f64, f64) -> f64, predicate: fn(f64, f64, f64) -> f64, direct: Direct) {
        self.param.push(Parameter::new_fun(ParamType::BoundaryCondition, value, predicate, direct));
    }
    pub fn add_concentrated_load_str(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.push(Parameter::new_str(ParamType::ConcentratedLoad, value, predicate, direct));
    }
    pub fn add_volume_load_str(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.push(Parameter::new_str(ParamType::VolumeLoad, value, predicate, direct));
    }
    pub fn add_volume_load_fun(&mut self, value: fn(f64, f64, f64) -> f64, predicate: fn(f64, f64, f64) -> f64, direct: Direct) {
        self.param.push(Parameter::new_fun(ParamType::VolumeLoad, value, predicate, direct));
    }
    pub fn add_surface_load_str(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.push(Parameter::new_str(ParamType::SurfaceLoad, value, predicate, direct));
    }
    pub fn add_surface_load_fun(&mut self, value: fn(f64, f64, f64) -> f64, predicate: fn(f64, f64, f64) -> f64, direct: Direct) {
        self.param.push(Parameter::new_fun(ParamType::SurfaceLoad, value, predicate, direct));
    }
    pub fn add_pressure_load_str(&mut self, value: &'a str, predicate: &'a str) {
        self.param.push(Parameter::new_str(ParamType::PressureLoad, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_pressure_load_fun(&mut self, value: fn(f64, f64, f64) -> f64, predicate: fn(f64, f64, f64) -> f64) {
        self.param.push(Parameter::new_fun(ParamType::PressureLoad, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_variable(&mut self, name: &'a str, value: f64) {
        self.variables.insert(name, value);
    }
}


pub trait FiniteElementMethod<'a>: Send + Sync {
    fn get_param(&self) -> Arc<&FEMParameter<'a>>;
    fn get_mesh(&self) -> Arc<&Mesh>;
    fn calc_fe_matrix(&self, i: usize) -> Result<Array2<f64>, FemError> {
        use crate::fem::fe::fe1d::FiniteElement1D;
        use crate::fem::fe::fe2d::FiniteElement2D;
        use crate::fem::fe::fe3d::FiniteElement3D;

        let x = self.get_mesh().get_fe_coord(i);
        let e = [self.get_param_value(&x, ParamType::YoungModulus)?, self.get_param_value(&x, ParamType::PoissonsRatio)?];
        let thk = self.get_param_value(&x, ParamType::Thickness)?;
        match self.get_mesh().fe_type {   
            FEType::FE1D2 => {
                let mut fe = fe::fe1d::FE1D2::new(e, thk, x);
                fe.generate()
            }
            FEType::FE2D3 | FEType::FE3D3S  => {
                let mut fe = fe::fe2d::FE2D3::new(e, thk, x, if self.get_mesh().fe_type == FEType::FE2D3 { false } else { true });
                fe.generate()
            }
            FEType::FE2D4 | FEType::FE3D4S => {
                let mut fe = fe::fe2d::FE2D4::new(e, thk, x, if self.get_mesh().fe_type == FEType::FE2D4 { false } else { true });
                fe.generate()
            }
            FEType::FE3D4 => {
                let mut fe = fe::fe3d::FE3D4::new(e, x);
                fe.generate()
            }
            FEType::FE3D8 => {
                let mut fe = fe::fe3d::FE3D8::new(e, x);
                fe.generate()
            }
        } 
    }
    fn set_global_matrix(&self, solver: &mut Mutex<impl Solver>) -> Result<(), FemError> {
        let msg = Mutex::new(Messenger::new("Generate global stiffness matrix", 1, self.get_mesh().num_fe as i64, 5));
        (0..self.get_mesh().num_fe).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
            msg.lock().unwrap().add_progress();
            let fe = self.calc_fe_matrix(i)?;
            for j in 0..fe.shape()[0] {
                for k in j..fe.shape()[1] {
                    let index1 = self.get_mesh().fe[[i, j / self.get_mesh().freedom]] * self.get_mesh().freedom + j % self.get_mesh().freedom;
                    let index2 = self.get_mesh().fe[[i, k / self.get_mesh().freedom]] * self.get_mesh().freedom + k % self.get_mesh().freedom;
                    solver.lock().unwrap().add_matrix_value(index1, index2, fe[[j, k]])?;
                    if j != k {
                        solver.lock().unwrap().add_matrix_value(index2, index1, fe[[j, k]])?;
                    }
                }
            }
            Ok(())
        })?;
        // msg.lock().unwrap().stop();
        Ok(())
    }
    fn num_results(&self) -> usize {
        match self.get_mesh().fe_type {
            FEType::FE1D2 => 3,
            FEType::FE2D3 | FEType::FE2D4 => 8,
            FEType::FE3D4 | FEType::FE3D8 => 15,
            FEType::FE3D3S | FEType::FE3D4S => 18,
        }
    }
    fn fun_names(&self) -> Vec<&str> {
        match self.get_mesh().fe_type {
            FEType::FE1D2 => vec![ "U", "Exx", "Sxx", "Si" ],
            FEType::FE2D3 | FEType::FE2D4 => vec![ "U", "V", "Exx", "Eyy", "Exy", "Sxx", "Syy", "Sxy", "Si" ],
            FEType::FE3D4 | FEType::FE3D8 => vec![ "U", "V", "W", "Exx", "Eyy", "Ezz", "Exy", "Exz", "Eyz", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz", "Si" ],
            FEType::FE3D3S | FEType::FE3D4S => vec![ "U", "V", "W", "Tx", "Ty", "Tz", "Exx", "Eyy", "Ezz", "Exy", "Exz", "Eyz", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz", "Si" ],
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
        let res = Mutex::new(Array2::zeros((self.num_results() + 1, self.get_mesh().num_vertex)));
        let counter = Mutex::new(Array1::<i32>::zeros(self.get_mesh().num_vertex));
        // Копирование результатов расчета (перемещений)
        for i in 0..self.get_mesh().freedom {
            for j in 0..self.get_mesh().num_vertex {
                res.lock().unwrap()[[i, j]] = u[j * self.get_mesh().freedom + i];
            }
        }
        // Вычисление деформаций, напряжений, ...
        let msg = Mutex::new(Messenger::new("Calculation of standard finite element results", 1, self.get_mesh().num_fe as i64, 5));
        (0..self.get_mesh().num_fe).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
            msg.lock().unwrap().add_progress();
            // Формируем вектор перемещений для текущего КЭ
            let mut fe_u = Array1::zeros(self.get_mesh().fe.shape()[1] * self.get_mesh().freedom);
            for j in 0..self.get_mesh().fe.shape()[1] {
                for k in 0..self.get_mesh().freedom {
                    fe_u[j * self.get_mesh().freedom + k] = u[self.get_mesh().freedom * self.get_mesh().fe[[i, j]] + k];
                }
            }
            let fe_res = self.calc_fe_results(i, fe_u)?; 
            for j in 0..self.num_results() - self.get_mesh().freedom {
                for k in 0..self.get_mesh().fe.shape()[1] {
                    res.lock().unwrap()[[j + self.get_mesh().freedom, self.get_mesh().fe[[i, k]]]] += fe_res[[j, k]];
                    if j == 0 {
                        counter.lock().unwrap()[self.get_mesh().fe[[i, k]]] += 1;
                    }
                }
            }
            Ok(())
        })?;
        let mut res = res.lock().unwrap();
        let counter = counter.lock().unwrap();
        // Осредняем результаты
        for i in self.get_mesh().freedom..self.num_results() {
            for j in 0..self.get_mesh().num_vertex {
                res[[i, j]] /= counter[j] as f64;
            }
        }
        // Вычисляем интенсивность напряжений
        let m_sqrt1_2 = 0.5 * 2_f64.sqrt();
        for i in 0..self.get_mesh().num_vertex {
            res[[self.num_results(), i]] = m_sqrt1_2 * match self.get_mesh().fe_type {
                FEType::FE1D2 => res[[2, i]].abs(),
                FEType::FE2D3 | FEType::FE2D4 => f64::sqrt(f64::powf(res[[5, i]] - res[[6, i]], 2.0) + 6.0 * (f64::powf(res[[7, i]], 2.0))),
                FEType::FE3D4 | FEType::FE3D8 => f64::sqrt(f64::powf(res[[9, i]] - res[[10, i]], 2.0) + f64::powf(res[[9, i]] - res[[11, i]], 2.0) +
                    f64::powf(res[[11, i]] - res[[12, i]], 2.0) + 6.0 * (f64::powf(res[[12, i]], 2.0) + f64::powf(res[[13, i]], 2.0) + f64::powf(res[[14, i]], 2.0))),
                FEType::FE3D3S | FEType::FE3D4S => f64::sqrt(f64::powf(res[[12, i]] - res[[13, i]], 2.0) + f64::powf(res[[12, i]] - res[[14, i]], 2.0) +
                    f64::powf(res[[13, i]] - res[[14, i]], 2.0) + 6.0 * (f64::powf(res[[15, i]], 2.0) + f64::powf(res[[16, i]], 2.0) + f64::powf(res[[17, i]], 2.0))),
            };
        }
        self.print_summary(&res);
        self.save_results(&res, res_name)
    }
    fn get_param_value(&self, x: &Array2<f64>, param_type: ParamType) -> Result<f64, FemError> {
        let mut val = 0.;
        for it in &self.get_param().param {
            if it.p_type == param_type {
                if !it.get_predicate(&x.row(0), &self.get_param().variables)? {
                    continue;
                }                
                val = it.get_value(&x.row(0), &self.get_param().variables)?;
                break;
            }
        }
        Ok(val)
    }
    fn calc_fe_results(&self, i: usize, u: Array1<f64>) -> Result<Array2<f64>, FemError> {
        use crate::fem::fe::fe1d::FiniteElement1D;
        use crate::fem::fe::fe2d::FiniteElement2D;
        use crate::fem::fe::fe3d::FiniteElement3D;

        let x = self.get_mesh().get_fe_coord(i);
        let e: [f64; 2] = [self.get_param_value(&x, ParamType::YoungModulus)?, self.get_param_value(&x, ParamType::PoissonsRatio)?];
        let thk = self.get_param_value(&x, ParamType::Thickness)?;
        match self.get_mesh().fe_type {   
            FEType::FE1D2 => {
                let fe = fe::fe1d::FE1D2::new(e, thk, x);
                fe.calc(&u)
            }
            FEType::FE2D3 | FEType::FE3D3S => {
                let fe = fe::fe2d::FE2D3::new(e, thk, x, if self.get_mesh().fe_type == FEType::FE2D3 { false } else { true });
                fe.calc(&u)
            }
            FEType::FE2D4 | FEType::FE3D4S => {
                let fe = fe::fe2d::FE2D4::new(e, thk, x, if self.get_mesh().fe_type == FEType::FE2D4 { false } else { true });
                fe.calc(&u)
            }
            FEType::FE3D4 => {
                let fe = fe::fe3d::FE3D4::new(e, x);
                fe.calc(&u)
            }
            FEType::FE3D8 => {
                let fe = fe::fe3d::FE3D8::new(e, x);
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
        if !write!(stream, "{}\n", self.get_mesh().fe_type).is_ok() {
            return Err(FemError::WriteFile);
        }
        if !write!(stream, "{}\n", self.get_mesh().num_vertex).is_ok() {
            return Err(FemError::WriteFile);
        }
        for i in 0..self.get_mesh().x.shape()[0] {
            for j in 0..self.get_mesh().x.shape()[1] {
                if !write!(stream, "{} ", self.get_mesh().x[[i, j]]).is_ok() {
                    return Err(FemError::WriteFile);
                }
            }
            if !write!(stream, "\n").is_ok() {
                return Err(FemError::WriteFile);
            }
        }
        if !write!(stream, "{}\n", self.get_mesh().num_fe).is_ok() {
            return Err(FemError::WriteFile);
        }
        for i in 0..self.get_mesh().fe.shape()[0] {
            for j in 0..self.get_mesh().fe.shape()[1] {
                if !write!(stream, "{} ", self.get_mesh().fe[[i, j]]).is_ok() {
                    return Err(FemError::WriteFile);
                }
            }
            if !write!(stream, "\n").is_ok() {
                return Err(FemError::WriteFile);
            }
        }
        if self.get_mesh().is_shell() {
            if !write!(stream, "0\n").is_ok() {
                return Err(FemError::WriteFile);
            }
        } else {
            if !write!(stream, "{}\n", self.get_mesh().num_be).is_ok() {
                return Err(FemError::WriteFile);
            }
            for i in 0..self.get_mesh().be.shape()[0] {
                for j in 0..self.get_mesh().be.shape()[1] {
                    if !write!(stream, "{} ", self.get_mesh().be[[i, j]]).is_ok() {
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
    fn volume_load_share(&self) -> Array1<f64> {
        match self.get_mesh().fe_type {   
            FEType::FE1D2 => array![ 0.5, 0.5 ],
            FEType::FE2D3 | FEType::FE3D3S => array![ 0.33333333333, 0.33333333333, 0.33333333333 ],
            FEType::FE2D4 | FEType::FE3D4S | FEType::FE3D4 => array![ 0.25, 0.25, 0.25, 0.25 ],
            FEType::FE3D8 => array![ 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125 ],
        } 
    }
    fn surface_load_share(&self) -> Array1<f64> {
        match self.get_mesh().fe_type {   
            FEType::FE1D2 => array![ 1.0 ],
            FEType::FE2D3 | FEType::FE2D4 => array![ 0.5, 0.5 ],
            FEType::FE3D3S | FEType::FE3D4 => array![ 0.333333333333, 0.333333333333, 0.333333333333 ],
            FEType::FE3D4S | FEType::FE3D8 => array![ 0.25, 0.25, 0.25, 0.25 ],
        } 
    }
    fn set_boundary_condition(&mut self, solver: &mut Mutex<impl Solver>) -> Result<(), FemError> {
        let msg = Mutex::new(Messenger::new("Using of boundary conditions", 1, self.get_mesh().num_vertex as i64, 5));
        (0..self.get_mesh().num_vertex).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
            msg.lock().unwrap().add_progress();
            for it in &self.get_param().param {
                if it.p_type == ParamType::BoundaryCondition {
                    let x = self.get_mesh().x.row(i);
                    if it.get_predicate(&x, &self.get_param().variables)? == false {
                        continue;
                    }
                    let val = it.get_value(&x, &self.get_param().variables)?;
                    if it.direct.contains(Direct::X) {
                        // solver.lock().unwrap().set_result_value((i + l) * self.mesh.freedom + 0, val)?;    
                        solver.lock().unwrap().set_result_value((i) * self.get_mesh().freedom + 0, val)?;    
                    }
                    if it.direct.contains(Direct::Y) && (self.get_mesh().is_2d() || self.get_mesh().is_3d()) {
                        solver.lock().unwrap().set_result_value((i) * self.get_mesh().freedom + 1, val)?;    
                    }
                    if it.direct.contains(Direct::Z) && self.get_mesh().is_3d() {
                        solver.lock().unwrap().set_result_value((i) * self.get_mesh().freedom + 2, val)?;    
                    }
                }
            }
            Ok(())    
        })?;
        // msg.lock().unwrap().stop();
        Ok(())
    }
    fn set_concentrated_load(&mut self, solver: &mut Mutex<impl Solver>) -> Result<(), FemError> {
        if self.get_param().find_parameter(ParamType::ConcentratedLoad) {
            let msg = Mutex::new(Messenger::new("Calculation of concentrated loads", 1, self.get_mesh().num_vertex as i64, 5));
            (0..self.get_mesh().num_vertex).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
                msg.lock().unwrap().add_progress();
                for it in &self.get_param().param {
                    if it.p_type == ParamType::ConcentratedLoad {
                        let x = self.get_mesh().x.row(i);
                        if it.get_predicate(&x, &self.get_param().variables)? == false {
                            continue;
                        }
                        let val = it.get_value(&x, &self.get_param().variables)?;
                        if it.direct.contains(Direct::X) {
                            solver.lock().unwrap().set_vector_value((i) * self.get_mesh().freedom + 0, val)?;    
                        }
                        if it.direct.contains(Direct::Y) && (self.get_mesh().is_2d() || self.get_mesh().is_3d()) {
                            solver.lock().unwrap().set_vector_value((i) * self.get_mesh().freedom + 1, val)?;    
                        }
                        if it.direct.contains(Direct::Z) && self.get_mesh().is_3d() {
                            solver.lock().unwrap().set_vector_value((i) * self.get_mesh().freedom + 2, val)?;    
                        }
                    }
                }
                Ok(())
            })?;
        }
        Ok(())
    }
    fn set_volume_load(&mut self, solver: &mut Mutex<impl Solver>) -> Result<(), FemError> {
        if self.get_param().find_parameter(ParamType::VolumeLoad) {
            let msg = Mutex::new(Messenger::new("Calculation of volume loads", 1, self.get_mesh().num_fe as i64, 5));
            (0..self.get_mesh().num_fe).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
                msg.lock().unwrap().add_progress();
                for it in &self.get_param().param {
                    if it.p_type == ParamType::VolumeLoad {
                        let x = self.get_mesh().get_fe_coord(i);
                        if !self.check_elem(&x, &it)? {
                            continue;
                        }                        
                        let share = self.volume_load_share() * self.get_mesh().fe_volume(i); 
                        for j in 0..self.get_mesh().fe.shape()[1] {
                            let val = it.get_value(&x.row(j), &self.get_param().variables)? * share[j];
                            if it.direct.contains(Direct::X) {
                                solver.lock().unwrap().add_vector_value(self.get_mesh().fe[[i, j]] * self.get_mesh().freedom + 0, val)?;    
                            }
                            if it.direct.contains(Direct::Y) && (self.get_mesh().is_2d() || self.get_mesh().is_3d()) {
                                solver.lock().unwrap().add_vector_value(self.get_mesh().fe[[i, j]] * self.get_mesh().freedom + 1, val)?;    
                            }
                            if it.direct.contains(Direct::Z) && self.get_mesh().is_3d() {
                                solver.lock().unwrap().add_vector_value(self.get_mesh().fe[[i, j]] * self.get_mesh().freedom + 2, val)?;    
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
        if self.get_param().find_parameter(ParamType::PressureLoad) || self.get_param().find_parameter(ParamType::SurfaceLoad) {
            let msg = Mutex::new(Messenger::new("Calculation of surface loads", 1, self.get_mesh().num_be as i64, 5));
            (0..self.get_mesh().num_be).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
                msg.lock().unwrap().add_progress();
                for it in &self.get_param().param {
                    if it.p_type == ParamType::PressureLoad || it.p_type == ParamType::SurfaceLoad {
                        let x = self.get_mesh().get_be_coord(i);
                        if !self.check_elem(&x, &it)? {
                            continue;
                        }                        
                        let share = self.surface_load_share() * self.get_mesh().be_volume(i); 
                        let normal = if it.p_type == ParamType::PressureLoad { self.get_mesh().be_normal(i) } else { array![1., 1., 1.] };
                        for j in 0..self.get_mesh().be.shape()[1] {
                            let val = it.get_value(&x.row(j), &self.get_param().variables)? * share[j];
                            if it.direct.contains(Direct::X) {
                                solver.lock().unwrap().add_vector_value(self.get_mesh().be[[i, j]] * self.get_mesh().freedom + 0, val * normal[0])?;    
                            }
                            if it.direct.contains(Direct::Y) && (self.get_mesh().is_2d() || self.get_mesh().is_3d()) {
                                solver.lock().unwrap().add_vector_value(self.get_mesh().be[[i, j]] * self.get_mesh().freedom + 1, val * normal[1])?;    
                            }
                            if it.direct.contains(Direct::Z) && self.get_mesh().is_3d() {
                                solver.lock().unwrap().add_vector_value(self.get_mesh().be[[i, j]] * self.get_mesh().freedom + 2, val * normal[2])?;    
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
            if param.get_predicate(&x.row(i), &self.get_param().variables)? == false {
                return Ok(false);
            }
        }
        Ok(true)
    }
    fn set_load(&mut self, solver: &mut Mutex<impl Solver>) -> Result<(), FemError> {
        self.set_concentrated_load(solver)?;
        self.set_volume_load(solver)?;
        self.set_surface_load(solver)
    }
}

pub struct FEM<'a> {
    mesh: &'a Mesh,
    param: &'a FEMParameter<'a>,
}

impl<'a> FiniteElementMethod<'a> for FEM<'a> {
    fn get_param(&self) -> Arc<&FEMParameter<'a>> {
        Arc::new(self.param)
    }
    fn get_mesh(&self) -> Arc<&Mesh> {
        Arc::new(self.mesh)
    }
}


impl<'a> FEM<'a> {
    pub fn new(mesh: &'a Mesh, param: &'a FEMParameter) -> Self {
        Self {mesh, param}
    }
    pub fn generate(&mut self, res_name: &str) -> Result<(), FemError> {
        rayon::ThreadPoolBuilder::new().num_threads(self.param.nthreads).build_global().unwrap();
        let time = Instant::now();
        let mut solver = Mutex::new(LzhSolver::new(&self.mesh));
        // let mut solver = Mutex::new(EnvSolver::new(&self.mesh));
        self.set_load(&mut solver)?;
        self.set_global_matrix(&mut solver)?;
        self.set_boundary_condition(&mut solver)?;
        self.calc_results(&solver.lock().unwrap().solve(self.param.eps)?, res_name)?;
        println!("Lead time: {:.2?}", time.elapsed());
        Ok(())
    }
 }

pub struct FEMPlasticity<'a> {
    mesh: &'a Mesh,
    param: &'a FEMParameter<'a>,
}

impl<'a> FiniteElementMethod<'a> for FEMPlasticity<'a> {
    fn get_param(&self) -> Arc<&FEMParameter<'a>> {
        Arc::new(self.param)
    }
    fn get_mesh(&self) -> Arc<&Mesh> {
        Arc::new(self.mesh)
    }
}


impl<'a> FEMPlasticity<'a> {
    pub fn new(mesh: &'a Mesh, param: &'a FEMParameter) -> Self {
        Self {mesh, param}
    }
    pub fn generate(&mut self, res_name: &str) -> Result<(), FemError> {
        rayon::ThreadPoolBuilder::new().num_threads(self.param.nthreads).build_global().unwrap();
        let time = Instant::now();
        let mut solver = Mutex::new(LzhSolver::new(&self.mesh));
        // let mut solver = Mutex::new(EnvSolver::new(&self.mesh));
        self.set_load(&mut solver)?;
        self.set_global_matrix(&mut solver)?;
        self.set_boundary_condition(&mut solver)?;
        self.calc_results(&solver.lock().unwrap().solve(self.param.eps)?, res_name)?;
        println!("Lead time: {:.2?}", time.elapsed());
        Ok(())
    }
 }
