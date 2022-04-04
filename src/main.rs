#[path = "core/fem.rs"]
mod fem;
#[path = "core/error.rs"]
mod error;

// 1. Изменить процедуру учета граничных условий (обрабатывать только ненулевые элементы)
// 2. Выборка индексов КЭ (переделать) (x)
// 3. Парсер - предкомпиляция (x)
// 4. Параллельное умножение матрицы на число (x)
// 5. unwrap() (Parser)
// 6. Sparse - попробовать крайт sprs
// 7. Запись в mesh-файл информации о связях (x)
// 8. Вывод на экран информации о погрешности при итерационном решении СЛАУ
// 9. Запись результатов (x)

use std::env;
use json;
use std::fs;
use fem::Direct;
use error::{ErrorCode, Error, error};

#[allow(dead_code)]
fn test_1d2() {
    let mesh_name = "D:/Work/python/pyfem/mesh/body1d.trpa";
    // let mesh_name = "/home/serg/work/python/pyfem/mesh/body1d.trpa";
    let mut fem: fem::FEM = match fem::FEM::new(mesh_name) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_num_threads(4);
    fem.set_young_modulus(203200.);
    // fem.set_poisons_ratio(0.33);
    fem.set_thickness(1.0);
    fem.set_eps(1.0e-6);
    fem.add_boundary_condition("0", "x == 0", Direct::X);
    fem.add_concentrated_load("1", "x == 10", Direct::X);
    match fem.generate("D:/Work/python/pyfem/mesh/body1d.res") {
        Err(err) => {
            println!("{}", err);
            return;
        }
        _ => println!("Done"),
    }
}

#[allow(dead_code)]
fn test_2d4() {
    let file_name = ("D:/Work/python/pyfem/mesh/console4.trpa", "D:/Work/python/pyfem/mesh/console4.res");
    // let file_name = ("/home/serg/work/python/pyfem/mesh/console4.trpa", "/home/serg/work/python/pyfem/mesh/console4.res");
    // let file_name = ("/home/serg/work/Qt/QFEM/mesh/console/console4.trpa", "/home/serg/work/Qt/QFEM/mesh/console/console4.res");
    let mut fem: fem::FEM = match fem::FEM::new(file_name.0) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_num_threads(4);
    fem.set_young_modulus(203200.0);
    fem.set_poisons_ratio(0.27);

    
    fem.set_thickness(1.0);
    fem.add_boundary_condition("0", "x == 0", Direct::X | Direct::Y);
    // fem.add_concentrated_load(String::from("-1"), String::from("x == 10"), Direct::Y);
    fem.add_volume_load("1", "", Direct::Y);
    match fem.generate(file_name.1) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        _ => println!("Done"),
    }
}

#[allow(dead_code)]
fn test_3d4() {
    let file_name = ("/home/serg/work/Qt/QFEM/mesh/balka.trpa", "/home/serg/work/Qt/QFEM/mesh/balka.res");
    // let file_name = ("D:/Work/Qt/QFEM/mesh/balka.trpa", "D:/Work/Qt/QFEM/mesh/balka.res");
    // let mesh_name = "D:/Work/Qt/QFEM/mesh/balka.trpa";
    let mut fem: fem::FEM = match fem::FEM::new(file_name.0) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_num_threads(4);
    fem.set_young_modulus(203200.);
    fem.set_poisons_ratio(0.27);
    fem.add_boundary_condition("0", "y == 0", Direct::X | Direct::Y | Direct::Z);
    // fem.add_concentrated_load(String::from("-1"), String::from("y == 4"), Direct::Y);
    fem.add_volume_load("-100", "", Direct::Y);
    match fem.generate(file_name.1) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        _ => println!("Done"),
    }
}

#[allow(dead_code)]
fn test_3d8() {
    let file_name = ("D:/work/python/pyfem/mesh/cube.trpa", "D:/work/python/pyfem/mesh/cube.res");
    // let file_name = ("/home/serg/work/python/pyfem/mesh/cube.trpa", "/home/serg/work/python/pyfem/mesh/cube.res");
    let mut fem: fem::FEM = match fem::FEM::new(file_name.0) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_num_threads(4);
    fem.set_young_modulus(203200.);
    fem.set_poisons_ratio(0.27);
    
    fem.add_boundary_condition("0", "z == 0", Direct::X | Direct::Y | Direct::Z);
    // fem.add_boundary_condition_fun(|_x, _y, _z| { 0. }, |_x, _y, z| { if z == 0. { true } else { false } }, Direct::X | Direct::Y | Direct::Z);
    
    // fem.add_concentrated_load("1", "z == 1", Direct::Z);
    // fem.add_volume_load("1", "", Direct::Z);
    // fem.add_surface_load("1", "z == 1", Direct::Z);
    
    // fem.add_pressure_load("1", "z == 1");

    // fem.add_surface_load_fun(|_x, _y, _z| { -0.5 }, |_x, _y, z| { if z == 1. { true } else { false } }, Direct::Z);
    
    // fem.add_volume_load_fun(|_x, _y, _z| { -0.5 }, |_x, _y, _z| { true }, Direct::Z);
    fem.add_volume_load("-0.5", "", Direct::Z);
    match fem.generate(file_name.1) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        _ => println!("Done"),
    }
}

#[allow(dead_code)]
fn test_shell_3() {
    let file_name = ("D:/work/python/pyfem/mesh/shell-tube-3.trpa", "D:/work/python/pyfem/mesh/shell-tube-3.res");
    // let file_name = ("/home/serg/work/python/pyfem/mesh/shell-tube-3.trpa", "/home/serg/work/python/pyfem/mesh/shell-tube-3.res");
    // let file_name = ("/home/serg/work/python/pyfem/mesh/shell4_1_0.trpa", "/home/serg/work/python/pyfem/mesh/shell4_1_0.res");
    // let file_name = ("D:/work/python/pyfem/mesh/shell4_1_0.trpa", "D:/work/python/pyfem/mesh/shell4_1_0.res");
    let mut fem: fem::FEM = match fem::FEM::new(file_name.0) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_num_threads(4);
    fem.set_young_modulus(203200.);
    fem.set_poisons_ratio(0.27);
    fem.set_thickness(0.0369);
    
    fem.add_boundary_condition_fun(|_x, _y, _z| { 0. }, |_x, _y, z| { if z == 0. || z == 4.014 { true } else { false } }, Direct::X | Direct::Y | Direct::Z);
    fem.add_pressure_load_fun(|_x, _y, _z| { 0.05 }, |_x, _y, _z| { true });

    match fem.generate(file_name.1) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        _ => println!("Done"),
    }
}


fn main() {
    let args: Vec<String> = env::args().collect();
    
    if args.len() < 2 {
        println!("Too few parameters!");
        return;
    }
    // test_1d2();
    // test_2d4();
    // test_3d4();
    // test_3d8();
    // test_shell_3();
    match read_json(args[1].as_str()) {
        Ok(_) => println!("Done"),
        Err(e) => println!("Error: {}", e),
    };
}

fn read_json(file_name: &str) -> Result<(), Error> {

    // let data = r#"
    // {
    //     "Mesh": "/home/serg/work/python/pyfem/mesh/cube.trpa",
    //     "Result": "/home/serg/work/python/pyfem/mesh/cube.res",
    //     "YoungModulus": 203200,
    //     "PoissonRatio": 0.27,
    //     "Threads": 4,
    //     "BoundaryConditions": [
    //         {
    //             "Value": "0", 
    //             "Predicate": "z == 0", 
    //             "Direct": "XYZ"
    //         }
    //     ],
    //     "VolumeLoad": [
    //         {
    //             "Value": "-0.5", 
    //             "Predicate": "", 
    //             "Direct": "Z"
    //         }
    //     ]
    // }"#;

    let data = match fs::read_to_string(file_name) {
        Ok(v) => v,
        Err(_) => return Err(error(ErrorCode::ReadFile)),
    };

    let v = match json::parse(data.as_str()) {
        Ok(v) => v,
        Err(e) => {
            println!("{}", e);
            return Err(error(ErrorCode::JsonError));
        }
    };

    let mesh_name = match v["Mesh"].as_str() {
        Some(v) => v,
        None => return Err(error(ErrorCode::MeshError)),
    };
    let mut fem: fem::FEM = fem::FEM::new(mesh_name)?;

    let nthreads = match v["Threads"].as_u32() {
        Some(v) => v,
        None => 1,
    };
    fem.set_num_threads(nthreads as usize);

    let thickness = match v["Thickness"].as_f64() {
        Some(v) => v,
        None => 1.0,
    };
    fem.set_thickness(thickness);

    let ym = match v["YoungModulus"].as_f64() {
        Some(v) => v,
        None => return Err(error(ErrorCode::YoungModulusError)),
    };
    fem.set_young_modulus(ym);

    let pr = match v["PoissonRatio"].as_f64() {
        Some(v) => v,
        None => return Err(error(ErrorCode::PoissonRatioError)),
    };
    fem.set_poisons_ratio(pr);

    for i in 0..v["BoundaryConditions"].len() {
        let direct = get_json_direct("BoundaryConditions", &v, i)?;
        let value: &str = match v["BoundaryConditions"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(error(ErrorCode::ValueError)),
        };
        let predicate: &str = match v["BoundaryConditions"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(error(ErrorCode::PredicateError)),
        };
        fem.add_boundary_condition(value, predicate, direct);
    }

    for i in 0..v["VolumeLoad"].len() {
        let direct = get_json_direct("VolumeLoad", &v, i)?;
        let value: &str = match v["VolumeLoad"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(error(ErrorCode::ValueError)),
        };
        let predicate: &str = match v["VolumeLoad"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(error(ErrorCode::PredicateError)),
        };
        fem.add_volume_load(value, predicate, direct);
    }

    for i in 0..v["ConcentratedLoad"].len() {
        let direct = get_json_direct("ConcentratedLoad", &v, i)?;
        let value: &str = match v["ConcentratedLoad"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(error(ErrorCode::ValueError)),
        };
        let predicate: &str = match v["ConcentratedLoad"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(error(ErrorCode::PredicateError)),
        };
        fem.add_concentrated_load(value, predicate, direct);
    }

    for i in 0..v["SurfaceLoad"].len() {
        let direct = get_json_direct("SurfaceLoad", &v, i)?;
        let value: &str = match v["SurfaceLoad"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(error(ErrorCode::ValueError)),
        };
        let predicate: &str = match v["SurfaceLoad"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(error(ErrorCode::PredicateError)),
        };
        fem.add_surface_load(value, predicate, direct);
    }

    for i in 0..v["PressureLoad"].len() {
        let value: &str = match v["PressureLoad"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(error(ErrorCode::ValueError)),
        };
        let predicate: &str = match v["PressureLoad"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(error(ErrorCode::PredicateError)),
        };
        fem.add_pressure_load(value, predicate);
    }

    let res_name = match v["Result"].as_str() {
        Some(v) => v,
        None => return Err(error(ErrorCode::ResultError)),
    };
    fem.generate(res_name)
}

fn get_json_direct<'a>(name: &'a str, v: &json::JsonValue, i: usize) -> Result<Direct, Error> {
    let mut direct: Option<Direct> = None;
    let direct_str = match v[name][i]["Direct"].as_str() {
        Some(v) => v,
        None => return Err(error(ErrorCode::DirectError)),    
    };
    if direct_str.contains("X") {
        direct = Some(Direct::X);
    }
    if direct_str.contains("Y") {
        direct = if direct.is_some() { Some(Direct::Y | direct.unwrap()) } else { Some(Direct::Y) };
    }
    if direct_str.contains("Z") {
        direct = if direct.is_some() { Some(Direct::Z | direct.unwrap()) } else { Some(Direct::Z) };
    }
    match direct {
        Some(v) => Ok(v),
        None => Err(error(ErrorCode::IncorrectDirectError)),
    }
}
