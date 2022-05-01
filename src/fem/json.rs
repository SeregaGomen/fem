use json;
use crate::fem;
use std::fs;
use crate::error::{ErrorCode, Error, error};
use crate::fem::Direct;


#[allow(dead_code)]
pub fn read_json(file_name: &str) -> Result<(), Error> {

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

    let eps = match v["Eps"].as_f64() {
        Some(v) => v,
        None => 1.0e-6,
    };
    fem.set_eps(eps);
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

#[allow(dead_code)]
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
