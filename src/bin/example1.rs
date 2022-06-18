extern crate fem;

use fem::fem::{generate, FEM};
use fem::fem::param::{Direct, FEMParameter};
use fem::fem::mesh::Mesh;
use fem::fem::error::FemError;

#[allow(dead_code)]
fn test_1d2(nthreads: usize) -> Result<(), FemError> {
    let file_name = ("data/body1d.mesh", "data/body1d.res");
    let mesh = Mesh::new(file_name.0)?;
    let mut param = FEMParameter::new();

    param.set_num_threads(nthreads);
    param.set_eps(1.0e-6);
    param.add_young_modulus_str("203200.", "");
    param.add_thickness_str("1.0", "");
    param.add_boundary_condition_str("0", "x == 0", Direct::X);
    param.add_concentrated_load_str("1", "x == 10", Direct::X);
    generate::<FEM>(&mesh, &param, file_name.1)
}

#[allow(dead_code)]
fn test_2d4(nthreads: usize) -> Result<(), FemError> {
    let file_name = ("data/console.mesh", "data/console.res");
    let mesh = Mesh::new(file_name.0)?;
    let mut param = FEMParameter::new();

    param.set_num_threads(nthreads);
    param.add_young_modulus_str("203200.0", "");
    param.add_poisons_ratio_str("0.27", "");
    param.add_thickness_str("1.0", "");
    param.add_boundary_condition_str("0", "x == 0", Direct::X | Direct::Y);
    param.add_concentrated_load_str("-1", "x == 10", Direct::Y);
    generate::<FEM>(&mesh, &param, file_name.1)
}

#[allow(dead_code)]
fn test_3d8(nthreads: usize) -> Result<(), FemError> {
    let file_name = ("data/cube.mesh", "data/cube.res");
    let mesh = Mesh::new(file_name.0)?;
    let mut param = FEMParameter::new();

    param.set_num_threads(nthreads);
    param.add_young_modulus_str("203200.0", "");
    param.add_poisons_ratio_str("0.27", "");
    
    // param.add_boundary_condition("0", "z == 0", Direct::X | Direct::Y | Direct::Z);
    param.add_boundary_condition_fun(|_x, _y, _z| { 0. }, |_x, _y, z| { if z == 0. { 1. } else { 0. } }, Direct::X | Direct::Y | Direct::Z);
    // param.add_volume_load_fun(|_x, _y, _z| { -0.5 }, |_x, _y, _z| { true }, Direct::Z);
    param.add_volume_load_str("-0.5", "", Direct::Z);
    generate::<FEM>(&mesh, &param, file_name.1)
}

#[allow(dead_code)]
fn test_shell_3(nthreads: usize) -> Result<(), FemError> {
    let file_name = ("data/shell-tube-3.mesh", "data/shell-tube-3.res");
    let mesh = Mesh::new(file_name.0)?;
    let mut param = FEMParameter::new();
    let fn_true = |_x: f64, _y: f64, _z: f64| { 1.0 };

    param.set_num_threads(nthreads);
    param.add_young_modulus_fun(|_x, _y, _z| { 6.5e+10 }, fn_true);
    param.add_poisons_ratio_fun(|_x, _y, _z| { 0.3 }, fn_true);
    param.add_thickness_fun(|_x, _y, _z| { 0.0045 }, fn_true);
    param.add_boundary_condition_fun(|_x, _y, _z| { 0. }, |_x, _y, z| { if z == 0. || z == 4.014 { 1. } else { 0. } }, Direct::X | Direct::Y | Direct::Z);
    param.add_pressure_load_fun(|_x, _y, _z| { 3.65E+07 }, fn_true);
    generate::<FEM>(&mesh, &param, file_name.1)
}


fn main() {
    // test_1d2(8);
    // test_2d4(8);
    // test_3d8(8);
    if let Err(e) = test_shell_3(8) {
        println!("\n\x1b[91mError: {}\x1b[0m", e);
    }
}


