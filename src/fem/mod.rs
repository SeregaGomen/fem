pub mod error;
pub mod mesh;
pub mod param;
mod fe;
mod util;
mod sparse;
mod solver;
mod parser;
mod msg;

use std::sync::{Arc, Mutex};
use rayon::prelude::*;
use ndarray::{Array1, Array2, prelude::*};
use std::time::Instant;
use error::FemError;
use mesh::Mesh;
use solver::{FemSolver, RussellSolver};
use fe::FEType;
use msg::Messenger;
use param::{ParamType, Parameter, FEMParameter, Direct};



pub trait FiniteElementMethod<'a>: Send + Sync {
    fn new(mesh: &'a Mesh, param: &'a FEMParameter) -> Self;
    fn generate(&mut self, res_name: &str) -> Result<(), FemError>;
    fn get_param(&self) -> Arc<&FEMParameter<'a>>;
    fn get_mesh(&self) -> Arc<&Mesh>;
    fn get_fe_param(&self, i: usize) -> Result<(Array2<f64>, [f64; 2], f64), FemError> {
        use std::fs::OpenOptions;
        use std::io::Write;        
        let mut file_e = OpenOptions::new().append(true).create(true).open("e.txt").expect("cannot open file");
        let mut file_thk = OpenOptions::new().append(true).create(true).open("thk.txt").expect("cannot open file");
        let e = match self.get_param_value(i, ParamType::YoungModulus)? {
            Some(val) => val,
            None => return Err(FemError::YoungModulusError),
        };
        let m = match self.get_param_value(i, ParamType::PoissonsRatio)? {
            Some(val) => val,
            None => {
                if self.get_mesh().is_1d() {
                    0.
                } else {
                    return Err(FemError::PoissonRatioError)
                }
            }
        };
        let thk = match self.get_param_value(i, ParamType::Thickness)? {
            Some(val) => val,
            None => {
                if self.get_mesh().is_3d() {
                    0.
                } else {
                    return Err(FemError::ThicknessError)
                }
            }
        };

        file_e.write_all(format!("{}\n", e).as_bytes()).expect("write failed");
        file_thk.write_all(format!("{}\n", thk).as_bytes()).expect("write failed");


        Ok((self.get_mesh().get_fe_coord(i), [e, m], thk))
    }
    fn calc_fe_matrix(&self, i: usize) -> Result<Array2<f64>, FemError> {
        use crate::fem::fe::fe1d::FiniteElement1D;
        use crate::fem::fe::fe2d::FiniteElement2D;
        use crate::fem::fe::fe3d::FiniteElement3D;

        let fe_param = self.get_fe_param(i)?;
        match self.get_mesh().fe_type {   
            FEType::FE1D2 => {
                let mut fe = fe::fe1d::FE1D2::new(fe_param.0, fe_param.1, fe_param.2);
                fe.generate()
            }
            FEType::FE2D3 | FEType::FE3D3S  => {
                let mut fe = fe::fe2d::FE2D3::new(fe_param.0, fe_param.1, fe_param.2, if self.get_mesh().fe_type == FEType::FE2D3 { false } else { true });
                fe.generate()
            }
            FEType::FE2D4 | FEType::FE3D4S => {
                let mut fe = fe::fe2d::FE2D4::new(fe_param.0, fe_param.1, fe_param.2, if self.get_mesh().fe_type == FEType::FE2D4 { false } else { true });
                fe.generate()
            }
            FEType::FE3D4 => {
                let mut fe = fe::fe3d::FE3D4::new(fe_param.0, fe_param.1);
                fe.generate()
            }
            FEType::FE3D8 => {
                let mut fe = fe::fe3d::FE3D8::new(fe_param.0, fe_param.1);
                fe.generate()
            }
        } 
    }
    fn set_global_matrix(&self, solver: &mut Mutex<impl FemSolver>) -> Result<(), FemError> {
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
            FEType::FE1D2 => 4,
            FEType::FE2D3 | FEType::FE2D4 => 9,
            FEType::FE3D4 | FEType::FE3D8 => 16,
            FEType::FE3D3S | FEType::FE3D4S => 19,
        }
        // self.fun_names().len()
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
    fn calc_results(&self, u: &Vec<f64>) -> Result<Array2<f64>, FemError> {
        let res = Mutex::new(Array2::<f64>::zeros((self.num_results() - self.get_mesh().freedom - 1, self.get_mesh().num_vertex)));
        let counter = Mutex::new(Array1::<i32>::zeros(self.get_mesh().num_vertex));
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
            for j in 0..self.num_results() - self.get_mesh().freedom - 1 {
                for k in 0..self.get_mesh().fe.shape()[1] {
                    res.lock().unwrap()[[j, self.get_mesh().fe[[i, k]]]] += fe_res[[j, k]];
                    if j == 0 {
                        counter.lock().unwrap()[self.get_mesh().fe[[i, k]]] += 1;
                    }
                }
            }
            Ok(())
        })?;

        let mut results = Array2::zeros((self.num_results(), self.get_mesh().num_vertex));
        let res = res.lock().unwrap();
        // Копирование результатов расчета (перемещений)
        for i in 0..self.get_mesh().freedom {
            for j in 0..self.get_mesh().num_vertex {
                results[[i, j]] = u[j * self.get_mesh().freedom + i];
            }
        }
        let counter = counter.lock().unwrap();
        // Осредняем результаты
        for i in 0..res.shape()[0] {
            for j in 0..res.shape()[1] {
                results[[i + self.get_mesh().freedom, j]] = res[[i, j]] / counter[j] as f64;
            }
        }
        // Вычисляем интенсивность напряжений
        let m_sqrt1_2 = 0.5 * 2_f64.sqrt();
        for i in 0..self.get_mesh().num_vertex {
            results[[self.num_results() - 1, i]] = m_sqrt1_2 * match self.get_mesh().fe_type {
                FEType::FE1D2 => results[[2, i]].abs(),
                FEType::FE2D3 | FEType::FE2D4 => f64::sqrt(f64::powf(results[[5, i]] - results[[6, i]], 2.0) + 6.0 * (f64::powf(results[[7, i]], 2.0))),
                FEType::FE3D4 | FEType::FE3D8 => f64::sqrt(f64::powf(results[[9, i]] - results[[10, i]], 2.0) + f64::powf(results[[9, i]] - results[[11, i]], 2.0) +
                    f64::powf(results[[11, i]] - results[[12, i]], 2.0) + 6.0 * (f64::powf(results[[12, i]], 2.0) + f64::powf(results[[13, i]], 2.0) + f64::powf(results[[14, i]], 2.0))),
                FEType::FE3D3S | FEType::FE3D4S => f64::sqrt(f64::powf(results[[12, i]] - results[[13, i]], 2.0) + f64::powf(results[[12, i]] - results[[14, i]], 2.0) +
                    f64::powf(results[[13, i]] - results[[14, i]], 2.0) + 6.0 * (f64::powf(results[[15, i]], 2.0) + f64::powf(results[[16, i]], 2.0) + f64::powf(results[[17, i]], 2.0))),
            };
        }
        Ok(results)
    }
    fn get_param_value(&self, i: usize, param_type: ParamType) -> Result<Option<f64>, FemError> {
        for it in &self.get_param().param {
            if it.p_type == param_type {
                if self.check_elem(&self.get_mesh().get_fe_coord(i), &it)? {
                    return Ok(Some(it.get_scalar_value(&self.get_mesh().get_fe_center(i), &self.get_param().variables)?));
                }
            }
        }
        Ok(None)
    }
    fn calc_fe_results(&self, i: usize, u: Array1<f64>) -> Result<Array2<f64>, FemError> {
        use crate::fem::fe::fe1d::FiniteElement1D;
        use crate::fem::fe::fe2d::FiniteElement2D;
        use crate::fem::fe::fe3d::FiniteElement3D;

        let fe_param = self.get_fe_param(i)?;
        match self.get_mesh().fe_type {   
            FEType::FE1D2 => {
                let fe = fe::fe1d::FE1D2::new(fe_param.0, fe_param.1, fe_param.2);
                fe.calc(&u)
            }
            FEType::FE2D3 | FEType::FE3D3S => {
                let fe = fe::fe2d::FE2D3::new(fe_param.0, fe_param.1, fe_param.2, if self.get_mesh().fe_type == FEType::FE2D3 { false } else { true });
                fe.calc(&u)
            }
            FEType::FE2D4 | FEType::FE3D4S => {
                let fe = fe::fe2d::FE2D4::new(fe_param.0, fe_param.1, fe_param.2, if self.get_mesh().fe_type == FEType::FE2D4 { false } else { true });
                fe.calc(&u)
            }
            FEType::FE3D4 => {
                let fe = fe::fe3d::FE3D4::new(fe_param.0, fe_param.1);
                fe.calc(&u)
            }
            FEType::FE3D8 => {
                let fe = fe::fe3d::FE3D8::new(fe_param.0, fe_param.1);
                fe.calc(&u)
            }
        } 
    }
    fn save_results(&self, res: &Array2<f64>, file_name: &str) -> Result<(), FemError> {
        use chrono::{Datelike, Timelike, Utc};
        use std::fs::File;
        use std::io::{BufWriter, prelude::*};
                
        let file = File::create(file_name)?;
        let mut stream = BufWriter::new(file);
        // Запись сигнатуры
        write!(stream, "FEM Solver Results File\n")?;
        // Вывод сетки
        write!(stream, "Mesh\n")?;
        write!(stream, "{}\n", self.get_mesh().fe_type)?;
        write!(stream, "{}\n", self.get_mesh().num_vertex)?;
        for i in 0..self.get_mesh().x.shape()[0] {
            for j in 0..self.get_mesh().x.shape()[1] {
                write!(stream, "{} ", self.get_mesh().x[[i, j]])?;
            }
            write!(stream, "\n")?;
        }
        write!(stream, "{}\n", self.get_mesh().num_fe)?;
        for i in 0..self.get_mesh().fe.shape()[0] {
            for j in 0..self.get_mesh().fe.shape()[1] {
                write!(stream, "{} ", self.get_mesh().fe[[i, j]])?;
            }
            write!(stream, "\n")?;
        }
        if self.get_mesh().is_shell() {
            write!(stream, "0\n")?;
        } else {
            write!(stream, "{}\n", self.get_mesh().num_be)?;
            for i in 0..self.get_mesh().be.shape()[0] {
                for j in 0..self.get_mesh().be.shape()[1] {
                    write!(stream, "{} ", self.get_mesh().be[[i, j]])?;
                }
                write!(stream, "\n")?;
            }
        }
        // Запись результатов расчета
        write!(stream, "Results\n")?;
        let now = Utc::now();
        write!(stream, "{:02}.{:02}.{:4} - {:02}:{:02}:{:02}\n",now.day(), now.month(), now.year(), now.hour(), now.minute(), now.second())?;
        write!(stream, "{}\n", self.num_results())?;
        for i in 0..res.shape()[0] {
            write!(stream, "{}\n", self.fun_names()[i])?;
            write!(stream, "0\n{}\n", res.shape()[1])?;
            for j in 0..res.shape()[1] {
                write!(stream, "{}\n", res[[i, j]])?;
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
    fn set_boundary_condition(&mut self, solver: &mut Mutex<impl FemSolver>) -> Result<(), FemError> {
        let msg = Mutex::new(Messenger::new("Using of boundary conditions", 1, self.get_mesh().num_vertex as i64, 5));
        (0..self.get_mesh().num_vertex).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
            msg.lock().unwrap().add_progress();
            for it in &self.get_param().param {
                if it.p_type == ParamType::BoundaryCondition {
                    if it.get_predicate(&self.get_mesh().x.row(i), &self.get_param().variables)? == true {
                        let val = it.get_scalar_value(&self.get_mesh().get_x_coord(i), &self.get_param().variables)?;
                        if it.direct.contains(Direct::X) {
                            // solver.lock().unwrap().set_result_value((i + l) * self.mesh.freedom + 0, val)?;    
                            solver.lock().unwrap().set_boundary_condition((i) * self.get_mesh().freedom + 0, val)?;    
                        }
                        if it.direct.contains(Direct::Y) && (self.get_mesh().is_2d() || self.get_mesh().is_3d()) {
                            solver.lock().unwrap().set_boundary_condition((i) * self.get_mesh().freedom + 1, val)?;    
                        }
                        if it.direct.contains(Direct::Z) && self.get_mesh().is_3d() {
                            solver.lock().unwrap().set_boundary_condition((i) * self.get_mesh().freedom + 2, val)?;    
                        }
                    }
                }
            }
            Ok(())    
        })?;
        // msg.lock().unwrap().stop();
        Ok(())
    }
    fn set_concentrated_load(&mut self, solver: &mut Mutex<impl FemSolver>) -> Result<(), FemError> {
        if self.get_param().find_parameter(ParamType::ConcentratedLoad) {
            let msg = Mutex::new(Messenger::new("Calculation of concentrated loads", 1, self.get_mesh().num_vertex as i64, 5));
            (0..self.get_mesh().num_vertex).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
                msg.lock().unwrap().add_progress();
                for it in &self.get_param().param {
                    if it.p_type == ParamType::ConcentratedLoad {
                        if it.get_predicate(&self.get_mesh().x.row(i), &self.get_param().variables)? == true {
                            let val = it.get_scalar_value(&self.get_mesh().get_x_coord(i), &self.get_param().variables)?;
                            if it.direct.contains(Direct::X) {
                                solver.lock().unwrap().add_vector_value((i) * self.get_mesh().freedom + 0, val)?;    
                            }
                            if it.direct.contains(Direct::Y) && (self.get_mesh().is_2d() || self.get_mesh().is_3d()) {
                                solver.lock().unwrap().add_vector_value((i) * self.get_mesh().freedom + 1, val)?;    
                            }
                            if it.direct.contains(Direct::Z) && self.get_mesh().is_3d() {
                                solver.lock().unwrap().add_vector_value((i) * self.get_mesh().freedom + 2, val)?;    
                            }
                        }
                    }
                }
                Ok(())
            })?;
        }
        Ok(())
    }
    fn set_volume_load(&mut self, solver: &mut Mutex<impl FemSolver>) -> Result<(), FemError> {
        if self.get_param().find_parameter(ParamType::VolumeLoad) {
            let msg = Mutex::new(Messenger::new("Calculation of volume loads", 1, self.get_mesh().num_fe as i64, 5));
            (0..self.get_mesh().num_fe).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
                msg.lock().unwrap().add_progress();
                for it in &self.get_param().param {
                    if it.p_type == ParamType::VolumeLoad {
                        if self.check_elem(&self.get_mesh().get_fe_coord(i), &it)? {
                            let share = self.volume_load_share() * self.get_mesh().fe_volume(i); 
                            let val = it.get_scalar_value(&self.get_mesh().get_fe_center(i), &self.get_param().variables)?;
                            for j in 0..self.get_mesh().fe.shape()[1] {
                                if it.direct.contains(Direct::X) {
                                    solver.lock().unwrap().add_vector_value(self.get_mesh().fe[[i, j]] * self.get_mesh().freedom + 0, val * share[j])?;    
                                }
                                if it.direct.contains(Direct::Y) && (self.get_mesh().is_2d() || self.get_mesh().is_3d()) {
                                    solver.lock().unwrap().add_vector_value(self.get_mesh().fe[[i, j]] * self.get_mesh().freedom + 1, val * share[j])?;    
                                }
                                if it.direct.contains(Direct::Z) && self.get_mesh().is_3d() {
                                    solver.lock().unwrap().add_vector_value(self.get_mesh().fe[[i, j]] * self.get_mesh().freedom + 2, val * share[j])?;    
                                }
                            }
                        }                        
                    }
                }
                Ok(())
            })?;
        }
        Ok(())
    }
    fn set_surface_load(&mut self, solver: &mut Mutex<impl FemSolver>) -> Result<(), FemError> {
        if self.get_param().find_parameter(ParamType::PressureLoad) || self.get_param().find_parameter(ParamType::SurfaceLoad) {
            let msg = Mutex::new(Messenger::new("Calculation of surface loads", 1, self.get_mesh().num_be as i64, 5));
            (0..self.get_mesh().num_be).into_par_iter().try_for_each(|i| -> Result<(), FemError> {
                msg.lock().unwrap().add_progress();
                for it in &self.get_param().param {
                    if it.p_type == ParamType::PressureLoad || it.p_type == ParamType::SurfaceLoad {
                        if self.check_elem(&self.get_mesh().get_be_coord(i), &it)? {
                            let share = self.surface_load_share() * self.get_mesh().be_volume(i); 
                            let normal = if it.p_type == ParamType::PressureLoad { self.get_mesh().be_normal(i) } else { array![1., 1., 1.] };
                            let val = it.get_scalar_value(&self.get_mesh().get_be_center(i), &self.get_param().variables)?;
                            for j in 0..self.get_mesh().be.shape()[1] {
                                if it.direct.contains(Direct::X) {
                                    solver.lock().unwrap().add_vector_value(self.get_mesh().be[[i, j]] * self.get_mesh().freedom + 0, val * normal[0] * share[j])?;    
                                }
                                if it.direct.contains(Direct::Y) && (self.get_mesh().is_2d() || self.get_mesh().is_3d()) {
                                    solver.lock().unwrap().add_vector_value(self.get_mesh().be[[i, j]] * self.get_mesh().freedom + 1, val * normal[1] * share[j])?;    
                                }
                                if it.direct.contains(Direct::Z) && self.get_mesh().is_3d() {
                                    solver.lock().unwrap().add_vector_value(self.get_mesh().be[[i, j]] * self.get_mesh().freedom + 2, val * normal[2] * share[j])?;    
                                }
                            }
                        }                        
                    }
                }
                Ok(())
            })?;
        }
        Ok(())
    }
    // Проверка предиката для всех узлов элемента
    // fn check_elem(&self, x: &Array2<f64>, param: &Parameter) -> Result<bool, FemError> {
    //     for i in 0..x.shape()[0] {
    //         if param.get_predicate(&x.row(i), &self.get_param().variables)? == false {
    //             return Ok(false);
    //         }
    //     }
    //     Ok(true)
    // }
    // Проверка предиката в центре элемента
    fn check_elem(&self, x: &Array2<f64>, param: &Parameter) -> Result<bool, FemError> {
        let mut cx = Array2::<f64>::zeros((1, 3));
        for j in 0..x.shape()[1] {
            let mut c = 0_f64;
            for i in 0..x.shape()[0] {
                c += x[[i, j]]
            }
            cx[[0, j]] = c / x.shape()[0] as f64;
        }
        if param.get_predicate(&cx.row(0), &self.get_param().variables)? == false {
            return Ok(false);
        }
        Ok(true)
    }
    fn set_load(&mut self, solver: &mut Mutex<impl FemSolver>) -> Result<(), FemError> {
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
    fn new(mesh: &'a Mesh, param: &'a FEMParameter) -> Self {
        Self {mesh, param}
    }
    fn get_param(&self) -> Arc<&FEMParameter<'a>> {
        Arc::new(self.param)
    }
    fn get_mesh(&self) -> Arc<&Mesh> {
        Arc::new(self.mesh)
    }
    fn generate(&mut self, res_name: &str) -> Result<(), FemError> {
        rayon::ThreadPoolBuilder::new().num_threads(self.param.nthreads).build_global().unwrap();
        let time = Instant::now();
        let mut solver = Mutex::new(RussellSolver::new(&self.mesh)?);
        self.set_boundary_condition(&mut solver)?;
        self.set_load(&mut solver)?;
        self.set_global_matrix(&mut solver)?;
        let res = self.calc_results(&solver.lock().unwrap().solve(self.param.eps)?)?;
        self.print_summary(&res);
        self.save_results(&res, res_name)?;
        println!("Lead time: {:.2?}", time.elapsed());
        Ok(())
    }
}


impl<'a> FEM<'a> {}

pub struct FEMPlasticity<'a> {
    mesh: &'a Mesh,
    param: &'a FEMParameter<'a>,
    iter_no: usize,
    res: Array2<f64>,
    ei: Mutex<(Vec<f64>, Vec<usize>)>,
    is_global_iteration_stop: Mutex<bool>,
    is_local_iteration_stop: Mutex<bool>,
}

impl<'a> FiniteElementMethod<'a> for FEMPlasticity<'a> {
    fn new(mesh: &'a Mesh, param: &'a FEMParameter) -> Self {
        Self {mesh, param, iter_no: 1, res: Array2::<f64>::zeros((0, 0)), ei: Mutex::new((Vec::new(), Vec::new())), is_global_iteration_stop: Mutex::new(false), is_local_iteration_stop: Mutex::new(false)}
    }
    fn get_param(&self) -> Arc<&FEMParameter<'a>> {
        Arc::new(self.param)
    }
    fn get_mesh(&self) -> Arc<&Mesh> {
        Arc::new(self.mesh)
    }
    fn generate(&mut self, res_name: &str) -> Result<(), FemError> {
        rayon::ThreadPoolBuilder::new().num_threads(self.param.nthreads).build_global().unwrap();
        let time = Instant::now();
        let mut solver = Mutex::new(RussellSolver::new(&self.mesh)?);

        // Учет краевых условий
        self.set_boundary_condition(&mut solver)?;
        // Предварительное вычисление компонент нагрузки
        self.set_load(&mut solver)?;
        
        let max_ssc = self.get_min_stress()?;
        let step = 0.05;
        let mut coef = 1.0;
        let mut add_count = 0.0;
        let mut count = 1;
        let mut is_loaded = false;

        *self.ei.lock().unwrap() = (vec![0.; self.get_mesh().fe.shape()[0]], vec![0; self.get_mesh().fe.shape()[0]]);
        // Итерационный процесс линеаризации упруго-пластической задачи
        while *self.is_global_iteration_stop.lock().unwrap() == false {
            *self.is_local_iteration_stop.lock().unwrap() = true;
            // Формирование ГМЖ
            self.set_global_matrix(&mut solver)?;
            let u = solver.lock().unwrap().solve(self.param.eps)?;
            if self.iter_no == 1 {
                // Вычисление интенсивности напряжений
                self.res = self.calc_results(&u)?;
                let max_si = util::get_max(self.res.row(self.res.shape()[0] - 1));
                if max_si > max_ssc {
                    // Задана слишком большая первоначальная нагрузка, уменьшаем ее на порядок
                    coef *= 0.1;
                    solver.lock().unwrap().mul_vector_value(0.1);
                    self.iter_no = 0;
                } else {
                    if is_loaded == false {
                        // Вычисляем поправочный коэффициент для "пропуска" упругой зоны
                        let load_factor = 0.95 * (max_ssc / max_si);
                        solver.lock().unwrap().mul_vector_value(load_factor);
                        coef *= load_factor;
                        self.iter_no = 0;
                        is_loaded = true;
                    } else {
                        // Устанавливаем нагрузку в значение "шаг по нагрузке"
                        solver.lock().unwrap().mul_vector_value(step);
                    }
                }
            } else {
                if *self.is_local_iteration_stop.lock().unwrap() {
                    self.res = &self.res + self.calc_results(&u)?;
                }
            }

            self.print_summary(&self.res);
            println!("Load: x {}", coef * (1.0 + add_count * step));
            println!("Iterate: {}\n", count);
            count += 1;
            self.iter_no += 1;
            if self.iter_no > 1 && *self.is_local_iteration_stop.lock().unwrap()  {
                add_count += 1.0;
                *self.is_local_iteration_stop.lock().unwrap() = false;
            }
            solver.lock().unwrap().clear_matrix();
        }
        self.save_results(&self.res, res_name)?;
        println!("Lead time: {:.2?}", time.elapsed());
        Ok(())
    }
    fn get_fe_param(&self, i: usize) -> Result<(Array2<f64>, [f64; 2], f64), FemError> {
        let mut e = match self.get_param_value(i, ParamType::YoungModulus)? {
            Some(val) => val,
            None => return Err(FemError::YoungModulusError),
        };
        let m = match self.get_param_value(i, ParamType::PoissonsRatio)? {
            Some(val) => val,
            None => {
                if self.get_mesh().is_1d() {
                    0.
                } else {
                    return Err(FemError::PoissonRatioError)
                }
            }
        };
        let thk = match self.get_param_value(i, ParamType::Thickness)? {
            Some(val) => val,
            None => {
                if self.get_mesh().is_3d() {
                    0.
                } else {
                    return Err(FemError::ThicknessError)
                }
            }
        };

        if self.iter_no > 1 {
            // Нелинейный случай
            let mut fe_si = self.res[[self.res.shape()[0] - 1, self.get_mesh().fe[[i, 0]]]];
            // Загружаем диаграмму деформирования, соответствующую текущему КЭ
            let mut ssc = &Vec::<[f64; 2]>::new();
            for it in &self.get_param().param {
                if it.p_type == ParamType::StressStrainCurve {
                    if !self.check_elem(&self.get_mesh().get_fe_coord(i), &it)? {
                        continue;
                    }                        
                    ssc = it.get_vector_value()?;
                    break;
                }
            }
            if ssc.len() == 0 {
                return Err(FemError::StressStrainCurveError);
            }
            // Определяем минимальную по КЭ интенсивность наряжений
            for j in 1..self.get_mesh().fe.shape()[1] {
                if fe_si < self.res[[self.res.shape()[0] - 1, self.get_mesh().fe[[i, j]]]] {
                    fe_si = self.res[[self.res.shape()[0] - 1, self.get_mesh().fe[[i, j]]]];            
                }
            }
            // Поиск в таблице свойств материала соответствующего напряжения
            let mut index = 0;
            if fe_si >= ssc[ssc.len() - 1][0] {
                // Достигнут предел текучести
                index = ssc.len() - 1;
                *self.is_global_iteration_stop.lock().unwrap() = true;
            } else {
                if fe_si >= ssc[1][0] {
                    for j in 1..ssc.len() {
                        if fe_si > ssc[j - 1][0] && fe_si <= ssc[j][0] { 
                            index = j;
                            break 
                        }
                    }
                }
            }
            let index_0 = self.ei.lock().unwrap().1[i];
            let new_e = if index != self.ei.lock().unwrap().1[i] {
                (ssc[index][0] / ssc[index][1] - ssc[index_0][0] / ssc[index_0][1]).abs()
            } else {
                if self.ei.lock().unwrap().0[i] == 0.0 { e } else { self.ei.lock().unwrap().0[i] }
            };
            e = new_e;
            // Проверка на изменение модуля упругости по сравнению с предыдущей итерацией
            if index != index_0 { *self.is_local_iteration_stop.lock().unwrap() = false }
            // Запоминаем рассчитанное значение модуля упругости и индекс
            self.ei.lock().unwrap().0[i] = new_e;
            self.ei.lock().unwrap().1[i] = index;
        }
        Ok((self.get_mesh().get_fe_coord(i), [e, m], thk))
    }
}

impl<'a> FEMPlasticity<'a> {
    fn get_min_stress(&self) -> Result<f64, FemError> {
        let mut res = f64::MAX;
        for i in &self.param.param {
            if i.p_type == ParamType::StressStrainCurve {
                let ssc = i.get_vector_value()?;
                if ssc.len() == 0 { return Err(FemError::StressStrainCurveError) }
                if ssc[1][0] < res { res = ssc[1][0] }
            }
        }
        Ok(res)
    }
}

pub fn generate<'a, F>(mesh: &'a Mesh, param: &'a FEMParameter, res_name: &str) -> Result<(), FemError> where F: FiniteElementMethod<'a> {
    let mut fem = F::new(mesh, param);
    fem.generate(res_name)
}