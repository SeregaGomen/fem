use ndarray::prelude::*;
use std::collections::HashMap;
use bitflags::bitflags;
use super::parser::Parser;
use super::error::FemError;


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
    Vector,
}

struct ParamValue<'a> {
    value_type: ParamValueType,
    value_str: &'a str,
    value_fun: Option<Box<dyn Fn(f64, f64, f64) -> f64 + Sync>>,
    value_vec: Vec<[f64; 2]>,
}

impl<'a> ParamValue<'a> {
    fn new_str(value: &'a str) -> Self {
        Self { value_type: ParamValueType::String, value_str: value, value_fun: None, value_vec: Vec::new() }
    }
    fn new_fun<F>(value: F) -> Self where F: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        Self { value_type: ParamValueType::Function, value_str: "", value_fun: Some(Box::new(value)), value_vec: Vec::new() }
    }
    fn new_vec(value: Vec<[f64; 2]>) -> Self {
        Self { value_type: ParamValueType::Vector, value_str: "", value_fun: None, value_vec: value }
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
    pub p_type: ParamType,
    pub direct: Direct,
    value: ParamValue<'a>,
    predicate: ParamValue<'a>,
}

impl<'a> Parameter<'a> {
    fn new_scalar_str(p_type: ParamType, val: &'a str, pred: &'a str, direct: Direct) -> Self {
        let value = ParamValue::new_str(val);
        let predicate = ParamValue::new_str(pred);
        Self { p_type, value, predicate, direct }
    }
    fn new_vector_str(p_type: ParamType, val: Vec<[f64; 2]>, pred: &'a str, direct: Direct) -> Self {
        let value = ParamValue::new_vec(val);
        let predicate = ParamValue::new_str(pred);
        Self { p_type, value, predicate, direct }
    }
    fn new_scalar_fun<V, P>(p_type: ParamType, val: V, pred: P, direct: Direct) -> Self where V: Fn(f64, f64, f64) -> f64 + Sync + 'static, P: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        let value = ParamValue::new_fun(val);
        let predicate = ParamValue::new_fun(pred);
        Self { p_type, value, predicate, direct }
    }
    fn new_vector_fun<F>(p_type: ParamType, val: Vec<[f64; 2]>, pred: F, direct: Direct) -> Self where F: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        let value = ParamValue::new_vec(val);
        let predicate = ParamValue::new_fun(pred);
        Self { p_type, value, predicate, direct }
    }
    pub fn get_scalar_value(&self, x: &Array1<f64>, variables: &HashMap<&'a str, f64>) -> Result<f64, FemError> {
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
                if let Some(fun) = &self.value.value_fun  {
                    Ok(fun(x[0], if x.len() > 1 { x[1] } else { 0. }, if x.len() > 2 { x[2] } else { 0. })) 
                } else {
                    Err(FemError::InternalError)
                }
            }
            _ => Err(FemError::IncorrectParamError)
        }
    }
    pub fn get_vector_value(&self) -> Result<&Vec<[f64; 2]>, FemError> {
        match self.value.value_type {
            ParamValueType::Vector => Ok(&self.value.value_vec),
            _ => Err(FemError::IncorrectParamError)
        }
    }
    pub fn get_predicate(&self, x: &ArrayView1<f64>, variables: &HashMap<&'a str, f64>) -> Result<bool, FemError> {
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
                let res = if let Some(fun) = &self.predicate.value_fun {
                    fun(x[0], if x.len() > 1 { x[1] } else { 0. }, if x.len() > 2 { x[2] } else { 0. } )
                } else {
                    return Err(FemError::InternalError)
                };
                Ok(if res == 1. {true} else {false})
            }
            _ => Err(FemError::IncorrectParamError)
        }
    }
}

pub struct FEMParameter<'a> {
    pub param: Vec<Parameter<'a>>,
    pub eps: f64,
    pub nthreads: usize,
    pub variables: HashMap<&'a str, f64>,
}

impl<'a> FEMParameter<'a> {
    pub fn new() -> Self {
        Self { param: Vec::new(), eps: 1.0e-6, nthreads: 1, variables: HashMap::new() }
    }
    pub fn find_parameter(&self, param: ParamType) -> bool {
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
    pub fn add_stress_strain_curve_str(&mut self, value: Vec<[f64; 2]>, predicate: &'a str) {
        self.param.push(Parameter::new_vector_str(ParamType::StressStrainCurve, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_stress_strain_curve_fun<P>(&mut self, value: Vec<[f64; 2]>, predicate: P) where P: Fn(f64, f64, f64) -> f64 + Sync + 'static, {
        self.param.push(Parameter::new_vector_fun(ParamType::StressStrainCurve, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_young_modulus_str(&mut self, value: &'a str, predicate: &'a str) {
        self.param.push(Parameter::new_scalar_str(ParamType::YoungModulus, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_young_modulus_fun<V, P>(&mut self, value: V, predicate: P) where V: Fn(f64, f64, f64) -> f64 + Sync + 'static, P: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        self.param.push(Parameter::new_scalar_fun(ParamType::YoungModulus, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_poisons_ratio_str(&mut self, value: &'a str, predicate: &'a str) {
        self.param.push(Parameter::new_scalar_str(ParamType::PoissonsRatio, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_poisons_ratio_fun<V, P>(&mut self, value: V, predicate: P) where V: Fn(f64, f64, f64) -> f64 + Sync + 'static, P: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        self.param.push(Parameter::new_scalar_fun(ParamType::PoissonsRatio, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_thickness_str(&mut self, value: &'a str, predicate: &'a str) {
        self.param.push(Parameter::new_scalar_str(ParamType::Thickness, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_thickness_fun<V, P>(&mut self, value: V, predicate: P) where V: Fn(f64, f64, f64) -> f64 + Sync + 'static, P: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        self.param.push(Parameter::new_scalar_fun(ParamType::Thickness, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn set_eps(&mut self, eps: f64) {
        self.eps = eps;
    }
    pub fn add_boundary_condition_str(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.push(Parameter::new_scalar_str(ParamType::BoundaryCondition, value, predicate, direct));
    }
    pub fn add_boundary_condition_fun<V, P>(&mut self, value: V, predicate: P, direct: Direct) where V: Fn(f64, f64, f64) -> f64 + Sync + 'static, P: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        self.param.push(Parameter::new_scalar_fun(ParamType::BoundaryCondition, value, predicate, direct));
    }
    pub fn add_concentrated_load_str(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.push(Parameter::new_scalar_str(ParamType::ConcentratedLoad, value, predicate, direct));
    }
    pub fn add_concentrated_load_fun<V, P>(&mut self, value: V, predicate: P, direct: Direct) where V: Fn(f64, f64, f64) -> f64 + Sync + 'static, P: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        self.param.push(Parameter::new_scalar_fun(ParamType::ConcentratedLoad, value, predicate, direct));
    }
    pub fn add_volume_load_str(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.push(Parameter::new_scalar_str(ParamType::VolumeLoad, value, predicate, direct));
    }
    pub fn add_volume_load_fun<V, P>(&mut self, value: V, predicate: P, direct: Direct) where V: Fn(f64, f64, f64) -> f64 + Sync + 'static, P: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        self.param.push(Parameter::new_scalar_fun(ParamType::VolumeLoad, value, predicate, direct));
    }
    pub fn add_surface_load_str(&mut self, value: &'a str, predicate: &'a str, direct: Direct) {
        self.param.push(Parameter::new_scalar_str(ParamType::SurfaceLoad, value, predicate, direct));
    }
    pub fn add_surface_load_fun<V, P>(&mut self, value: V, predicate: P, direct: Direct) where V: Fn(f64, f64, f64) -> f64 + Sync + 'static, P: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        self.param.push(Parameter::new_scalar_fun(ParamType::SurfaceLoad, value, predicate, direct));
    }
    pub fn add_pressure_load_str(&mut self, value: &'a str, predicate: &'a str) {
        self.param.push(Parameter::new_scalar_str(ParamType::PressureLoad, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_pressure_load_fun<V, P>(&mut self, value: V, predicate: P) where V: Fn(f64, f64, f64) -> f64 + Sync + 'static, P: Fn(f64, f64, f64) -> f64 + Sync + 'static {
        self.param.push(Parameter::new_scalar_fun(ParamType::PressureLoad, value, predicate, Direct::X | Direct::Y | Direct::Z));
    }
    pub fn add_variable(&mut self, name: &'a str, value: f64) {
        self.variables.insert(name, value);
    }
}
