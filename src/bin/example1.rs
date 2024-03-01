extern crate fem;

use fem::fem::{generate, FEM, FEMPlasticity};
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

#[allow(dead_code)]
fn test_tank(nthreads: usize) -> Result<(), FemError> {
    // let file_name = ("data/tank.mesh", "data/tank.res");
    let file_name = ("data/tank_1_4.mesh", "data/tank_1_4.res");
    let mesh = Mesh::new(file_name.0)?;
    let mut param = FEMParameter::new();
    let fn_true = |_x: f64, _y: f64, _z: f64| { 1.0 };
    let eps = 0.01;
    let e1 = 6.5e+10;
    let e2 = 7.3e+10;
    let c = 1.454;
    let cx_bot = 20.7657;
    let cx_top = -8.5497;
    let fi_b = -0.872665_f64;
    let fi_t = -2.26893_f64;
    let h = 0.06;
    let k2_bot = 0.0520196;
    let k2_top = 0.0520196;
    let l = 12.216;
    let l3 = 1.654;
    let l4 = 1.09;
    let p = 142196.0;
    let r = 2.5;

    param.set_num_threads(nthreads);
    param.add_young_modulus_fun(move |_x, _y, _z| { e1 }, move|x, y, z| { 
        if ((r * r - ((x - c) * (x - c)  + y * y + z * z)).abs() <= eps && x <= (r * fi_t.cos() + c)) || ((r * r - ((x - l + c) * (x - l + c) + y * y + z * z)).abs() <= eps && x >= (r * fi_b.cos() + l - c)) { 1.0 } else { 0.0 }
    });
    param.add_young_modulus_fun(move |_x, _y, _z| { e2 }, move|_, _, _| { 
        1.0
    });
    param.add_poisons_ratio_fun(|_x, _y, _z| { 0.3 }, fn_true);

    param.add_thickness_fun(|_x, _y, _z| { 0.0046 }, move |x, y, z| {
        if (((r * r - ((x - c) * (x - c) + y * y + z * z)).abs() <= eps) && (x <= (r * fi_t.cos() + c))) || (((r * r - ((x - l + c) * (x - l + c) + y * y + z * z)).abs() <= eps) && (x >= (r * fi_b.cos()) + l - c)) { 1.0 } else { 0.0 }
    });
    param.add_thickness_fun(|_x, _y, _z| { 0.05 }, move |x, _y, _z| {
        if ((x >= (r * fi_t.cos() + c)) && (x <= 0.)) || ((x >= l) && (x <= (r * fi_b.cos() + l - c))) || ((x >= 4. * l3 - h / 2.) && (x <= 4. * l3 + h / 2.)) { 1.0 } else { 0.0 }
    });
    param.add_thickness_fun(|_x, _y, _z| { 0.0255 }, move |x, _y, _z| {
        if ((x >= l3 - h / 2.0) && (x <= l3 + h / 2.)) || ((x >= 2. * l3 - h / 2.) && (x <= 2. * l3 + h / 2.)) || ((x >= 5. * l3 - h / 2.) && (x <= 5. * l3 + h / 2.)) || ((x >= 6. * l3 - h / 2.) && (x <= 6. * l3 + h / 2.)) || ((x >= 6. * l3 - h / 2. + l4) && (x <= 6. * l3 + h / 2. + l4)) { 1.0 } else { 0.0 }
    });
    param.add_thickness_fun(|_x, _y, _z| { 0.04 }, move |x, _y, _z| {
        if (x >= 3. * l3 - h) && (x <= 3. * l3 + h) { 1.0 } else { 0.0 }
    });
    param.add_thickness_fun(|_x, _y, _z| { 0.0045 }, move |x, _y, _z| {
        if (x >= 0. && x <= (l3 - h / 2.)) || (x >= (l3 + h / 2.) && x <= (2. * l3 - h / 2.)) || (x >= (2. * l3 + h / 2.) && x <= (3. * l3 - h)) || (x >= (4. * l3 + h / 2.) && x <= (5. * l3 - h / 2.)) || (x >= (5. * l3 + h / 2.) && x <= (6. * l3 - h / 2.)) { 1.0 } else { 0.0 }
    });
    param.add_thickness_fun(|_x, _y, _z| { 0.0046 }, move |x, _y, _z| {
        if (x >= (3. * l3 + h)) && (x <= (4. * l3 - h / 2.)) { 1.0 } else { 0.0 }
    });
    param.add_thickness_fun(|_x, _y, _z| { 0.0052 }, move |x, _y, _z| {
        if (x >= (6. * l3 + h / 2.) && x <= (6. * l3 - h / 2. + l4)) || (x >= (6. * l3 + h / 2. + l4) && x <= l) { 1.0 } else { 0.0 }
    });
    param.add_thickness_fun(|_x, _y, _z| { 0.0143 }, move |x, _y, _z| {
        if x < 0. { 1.0 } else { 0.0 }
    });
    param.add_thickness_fun(|_x, _y, _z| { 0.016 }, fn_true);
    
    param.add_boundary_condition_fun(|_, _, _| { 0. }, move |x, _, _| { if (x - 14.338).abs() < eps { 1.0 } else { 0.0 } }, Direct::X | Direct::Y | Direct::Z);
    param.add_boundary_condition_fun(|_, _, _| { 0. }, move |_, y, _| { if y.abs() <= eps { 1.0 } else { 0.0 } }, Direct::Y);
    param.add_boundary_condition_fun(|_, _, _| { 0. }, move |_, _, z| { if z.abs() <= eps { 1.0 } else { 0.0 } }, Direct::Z);

    param.add_pressure_load_fun(move |_x, _y, _z| { p }, move |x, _, _| { if x >= 0. && x <= l { 1.0 } else { 0.0 } });
    param.add_pressure_load_fun(move |_x, _y, _z| { p }, move |x, y, z| { if (r * r - ((x - c) * (x - c) + y * y + z * z)).abs() <= eps && x <= (r * fi_t.cos() + c) { 1.0 } else { 0.0 } });
    param.add_pressure_load_fun(move |_x, _y, _z| { p }, move |x, y, z| { if (r * r - ((x - l + c) * (x - l + c) + y * y + z * z)).abs() <= eps && x >= (r * fi_b.cos() + l - c) { 1.0 } else { 0.0 } });
    param.add_pressure_load_fun(move |_x, _y, _z| { p }, move |x, y, z| { if (x >= (r * fi_t.cos() + c) && x <= 0.) && (y * y + z * z - k2_top * (x - cx_top) * (x - cx_top)).abs() < eps { 1.0 } else { 0.0 } });
    param.add_pressure_load_fun(move |_x, _y, _z| { p }, move |x, y, z| { if (x >= l && x <= (r * fi_b.cos() + l - c)) && (y * y + z * z  - k2_bot * (x - cx_bot) * (x - cx_bot)).abs() < eps { 1.0 } else { 0.0 } });

    param.add_stress_strain_curve_fun(vec![[1.27E+08, 0.002], [1.57E+08, 0.004], [1.77E+08, 0.008], [1.96E+08, 0.02], [3.14E+08, 0.12]], move |x, y, z| { 
        if ((r * r - ((x - c) * (x - c)  + y * y + z * z)).abs() <= eps && x <= (r * fi_t.cos() + c)) || ((r * r - ((x - l + c) * (x - l + c) + y * y + z * z)).abs() <= eps && x >= (r * fi_b.cos() + l - c)) { 1.0 } else { 0.0 } 
    });
    param.add_stress_strain_curve_fun(vec![[1.96E+08, 0.003], [2.55E+08, 0.005], [2.75E+08, 0.006], [3.14E+08, 0.015], [3.92E+08, 0.0725]], fn_true);

    generate::<FEMPlasticity>(&mesh, &param, file_name.1)
    // generate::<FEM>(&mesh, &param, file_name.1)
}


fn main() {
    // if let Err(e) = test_shell_3(7) {
    //     println!("\n\x1b[91m{}\x1b[0m", e);
    // }
    if let Err(e) = test_tank(7) {
        println!("\n\x1b[91m{}\x1b[0m", e);
    }
    // if let Err(e) = test_1d2(1) {
    //     println!("\n\x1b[91m{}\x1b[0m", e);
    // }
    // if let Err(e) = test_3d8(1) {
    //     println!("\n\x1b[91m{}\x1b[0m", e);
    // }
}


