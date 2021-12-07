mod error;
mod fe;
mod util;
mod mesh;
mod sparse;
mod solver;
mod parser;
mod msg;

use bitflags::bitflags;
use ndarray::{Array1, Array2, prelude::*};
use std::time::Instant;
use crate::error::Error;
use mesh::Mesh;
use solver::Solver;
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
    fn get_value(&self, x: &ArrayView1<f64>) -> Result<f64, Error> {
        match self.value_fun {
            None => {
                let var_name = array!["x", "y", "z"];
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
    fn get_predicate(&self, x: &ArrayView1<f64>) -> Result<bool, Error> {
        match self.predicate_fun {
            None => {
                if self.predicate.len() == 0 {
                    Ok(true)
                }
                else {
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
}

impl<'a> FEMParameter<'a> {
    fn new() -> Self {
        Self { param: Vec::new(), e: [0., 0.], thk: 0., eps: 1.0e-6 }
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
    pub fn new(mesh_name: &str) -> Result<Self, Error> {
        let mesh = Mesh::new(mesh_name)?;
        Ok(Self {mesh, param: FEMParameter::new()})
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
    pub fn generate(&mut self) -> Result<(), Error> {
        let time = Instant::now();
        let mut solver = Solver::new(&self.mesh);
        let mut msg = Messenger::new(String::from("Generate global stiffness matrix"), 1, self.mesh.num_fe as i64, 5);
        for i in 0..self.mesh.num_fe {
            msg.add_progress();
            //print!("\rGenerate global stiffness matrix: {}", i);
            // let mut fe = FE::new(self.mesh.fe_type, self.param.e, self.param.thk, self.mesh.ge_fe_coord(i));
            // let lm = fe.generate()?;
            let lm = self.calc_fe_matrix(i)?;
            //println!("{:?}", lm);
            solver.add_local_matrix(&lm, &self.mesh.fe, i)?;
        }
        self.set_concentrated_load(&mut solver)?;
        self.set_volume_load(&mut solver)?;
        self.set_surface_load(&mut solver)?;
        //self.set_pressure_load(&mut solver)?;
        self.set_boundary_condition(&mut solver)?;
        // let res = self.calc_results(&solver.cg_solve(self.param.eps)?);
        // let res = self.calc_results(&solver.cho_solve(self.param.eps)?);
        self.print_summary(&self.calc_results(&solver.lzh_solve(self.param.eps)?)?);
        println!("Lead time: {:.2?}", time.elapsed());
        Ok(())
    }
    fn num_results(&self) -> usize {
        match self.mesh.fe_type {
            FEType::FE1D2 => 3,
            FEType::FE2D3 | FEType::FE2D4 => 8,
            FEType::FE3D4 | FEType::FE3D8 => 15,
        }
    }
    fn fun_names(&self) -> Vec<&str> {
        match self.mesh.fe_type {
            FEType::FE1D2 => vec![ "U", "Exx", "Sxx" ],
            FEType::FE2D3 | FEType::FE2D4 => vec![ "U", "V", "Exx", "Eyy", "Exy", "Sxx", "Syy", "Sxy" ],
            FEType::FE3D4 | FEType::FE3D8 => vec![ "U", "V", "W", "Exx", "Eyy", "Ezz", "Exy", "Exz", "Eyz", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz" ],
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
    fn calc_results(&self, u: &Array1<f64>) -> Result<Array2<f64>, Error> {
        let fe_size = self.mesh.fe.shape()[1];
        let freedom = self.mesh.freedom;
        let mut res: Array2<f64> = Array2::zeros((self.num_results(), self.mesh.num_vertex));
        let mut counter: Array1<i16> = Array1::zeros(self.mesh.num_vertex);
        // Копирование результатов расчета (перемещений)
        for i in 0..freedom {
            for j in 0..self.mesh.num_vertex {
                res[[i, j]] = u[[j * freedom + i]];
            }
        }
        // Вычисление деформаций, напряжений, ...
        let mut msg = Messenger::new(String::from("Calculation of standard finite element results"), 1, self.mesh.num_fe as i64, 5);
        for i in 0..self.mesh.num_fe {
            msg.add_progress();
            // Формируем вектор перемещений для текущего КЭ
            let mut fe_u = Array1::zeros(fe_size * freedom);
            for j in 0..fe_size {
                for k in 0..freedom {
                    fe_u[j * freedom + k] = u[freedom * self.mesh.fe[[i, j]] + k];
                }
            }
            // println!("\n{:?}", fe_u);
            let fe_res = self.calc_fe_results(i, fe_u)?; 
            // println!("\n{:?}", fe_res);
            for j in 0..self.num_results() - freedom {
                for k in 0..fe_size {
                    res[[j + freedom, self.mesh.fe[[i, k]]]] += fe_res[[j, k]];
                    if j == 0 {
                        counter[self.mesh.fe[[i, k]]] += 1;
                    }
                }
            }
        }
        // Осредняем результаты
        for i in freedom..self.num_results() {
            for j in 0..self.mesh.num_vertex {
                res[[i, j]] /= counter[j] as f64;
            }
        }
        Ok(res)
    }
    fn set_boundary_condition(&mut self, solver: &mut Solver) -> Result<(), Error> {
        let mut msg = Messenger::new(String::from("Using of boundary conditions"), 1, self.mesh.num_vertex as i64, 5);
        for i in 0..self.mesh.num_vertex {
            msg.add_progress();
            for it in &self.param.param {
                if it.p_type == ParamType::BoundaryCondition {
                    let x = self.mesh.x.row(i);
                    if it.get_predicate(&x)? == false {
                        continue;
                    }
                    let val = it.get_value(&x)?;
                    if it.direct.contains(Direct::X) {
                        solver.set_boundary_condition(i * self.mesh.freedom + 0, val)?;    
                    }
                    if it.direct.contains(Direct::Y) && (self.mesh.is_2d() || self.mesh.is_3d()) {
                        solver.set_boundary_condition(i * self.mesh.freedom + 1, val)?;    
                    }
                    if it.direct.contains(Direct::Z) && self.mesh.is_3d() {
                        solver.set_boundary_condition(i * self.mesh.freedom + 2, val)?;    
                    }
                }
            }
        }
        Ok(())
    }
    fn set_concentrated_load(&mut self, solver: &mut Solver) -> Result<(), Error> {
        if self.param.find_parameter(ParamType::ConcentratedLoad) {
            let mut msg = Messenger::new(String::from("Calculation of concentrated loads"), 1, self.mesh.num_vertex as i64, 5);
            for i in 0..self.mesh.num_vertex {
                msg.add_progress();
                for it in &self.param.param {
                    if it.p_type == ParamType::ConcentratedLoad {
                        let x = self.mesh.x.row(i);
                        if it.get_predicate(&x)? == false {
                            continue;
                        }
                        let val = it.get_value(&x)?;
                        if it.direct.contains(Direct::X) {
                            solver.set_load(i * self.mesh.freedom + 0, val)?;    
                        }
                        if it.direct.contains(Direct::Y) && (self.mesh.is_2d() || self.mesh.is_3d()) {
                            solver.set_load(i * self.mesh.freedom + 1, val)?;    
                        }
                        if it.direct.contains(Direct::Z) && self.mesh.is_3d() {
                            solver.set_load(i * self.mesh.freedom + 2, val)?;    
                        }
                    }
                }
            }
        }
        Ok(())
    }
    fn set_volume_load(&mut self, solver: &mut Solver) -> Result<(), Error> {
        if self.param.find_parameter(ParamType::VolumeLoad) {
            let mut msg = Messenger::new(String::from("Calculation of volume loads"), 1, self.mesh.num_fe as i64, 5);
            for i in 0..self.mesh.num_fe {
                msg.add_progress();
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
                                solver.add_load(self.mesh.fe[[i, j]] * self.mesh.freedom + 0, val)?;    
                            }
                            if it.direct.contains(Direct::Y) && (self.mesh.is_2d() || self.mesh.is_3d()) {
                                solver.add_load(self.mesh.fe[[i, j]] * self.mesh.freedom + 1, val)?;    
                            }
                            if it.direct.contains(Direct::Z) && self.mesh.is_3d() {
                                solver.add_load(self.mesh.fe[[i, j]] * self.mesh.freedom + 2, val)?;    
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
    fn set_surface_load(&mut self, solver: &mut Solver) -> Result<(), Error> {
        if self.param.find_parameter(ParamType::PressureLoad) || self.param.find_parameter(ParamType::SurfaceLoad) {
            let mut msg = Messenger::new(String::from("Calculation of surface loads"), 1, self.mesh.num_be as i64, 5);
            for i in 0..self.mesh.num_be {
                msg.add_progress();
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
                                solver.add_load(self.mesh.be[[i, j]] * self.mesh.freedom + 0, val * normal[0])?;    
                            }
                            if it.direct.contains(Direct::Y) && (self.mesh.is_2d() || self.mesh.is_3d()) {
                                solver.add_load(self.mesh.be[[i, j]] * self.mesh.freedom + 1, val * normal[1])?;    
                            }
                            if it.direct.contains(Direct::Z) && self.mesh.is_3d() {
                                solver.add_load(self.mesh.be[[i, j]] * self.mesh.freedom + 2, val * normal[2])?;    
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
    fn check_elem(&self, x: &Array2<f64>, param: &Parameter) -> Result<bool, Error> {
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
            FEType::FE2D3 => array![ 0.33333333333, 0.33333333333, 0.33333333333 ],
            FEType::FE2D4 | FEType::FE3D4 => array![ 0.25, 0.25, 0.25, 0.25 ],
            FEType::FE3D8 => array![ 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125 ],
        } 
    }
    fn surface_load_share(&self) -> Array1<f64> {
        match self.mesh.fe_type {   
            FEType::FE1D2 => array![ 1.0 ],
            FEType::FE2D3 | FEType::FE2D4 => array![ 0.5, 0.5 ],
            FEType::FE3D4 => array![ 0.333333333333, 0.333333333333, 0.333333333333 ],
            FEType::FE3D8 => array![ 0.25, 0.25, 0.25, 0.25 ],
        } 
    }
    fn calc_fe_matrix(&self, i: usize) -> Result<Array2<f64>, Error> {
        use crate::fem::fe::fe1d::FiniteElement1D;
        use crate::fem::fe::fe2d::FiniteElement2D;
        use crate::fem::fe::fe3d::FiniteElement3D;

        match self.mesh.fe_type {   
            FEType::FE1D2 => {
                let mut fe = fe::fe1d::FE1D2::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i));
                fe.generate()
            }
            FEType::FE2D3 => {
                let mut fe = fe::fe2d::FE2D3::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i));
                fe.generate()
            }
            FEType::FE2D4 => {
                let mut fe = fe::fe2d::FE2D4::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i));
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
    fn calc_fe_results(&self, i: usize, u: Array1<f64>) -> Result<Array2<f64>, Error> {
        use crate::fem::fe::fe1d::FiniteElement1D;
        use crate::fem::fe::fe2d::FiniteElement2D;
        use crate::fem::fe::fe3d::FiniteElement3D;

        match self.mesh.fe_type {   
            FEType::FE1D2 => {
                let fe = fe::fe1d::FE1D2::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i));
                fe.calc(&u)
            }
            FEType::FE2D3 => {
                let fe = fe::fe2d::FE2D3::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i));
                fe.calc(&u)
            }
            FEType::FE2D4 => {
                let fe = fe::fe2d::FE2D4::new(self.param.e, self.param.thk, self.mesh.get_fe_coord(i));
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
}