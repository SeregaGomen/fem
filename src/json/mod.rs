use json;
use std::fs;
use crate::fem::{FEM, FEMPlasticity,  error::FemError, generate};
use crate::fem::mesh::Mesh;
use crate::fem::param::{Direct, FEMParameter};

#[allow(dead_code)]
pub fn read_json(file_name: &str) -> Result<(), FemError> {
    let mut is_plasticity = false;
    let data = fs::read_to_string(file_name)?;
    let v = json::parse(data.as_str())?;
    let mesh_name = match v["Mesh"].as_str() {
        Some(v) => v,
        None => return Err(FemError::MeshError),
    };
    let mesh = Mesh::new(mesh_name)?;
    let mut param = FEMParameter::new();

    let eps = match v["Eps"].as_f64() {
        Some(v) => v,
        None => 1.0e-6,
    };
    param.set_eps(eps);
    let nthreads = match v["Threads"].as_u32() {
        Some(v) => v,
        None => 1,
    };
    param.set_num_threads(nthreads as usize);

    for i in 0..v["Thickness"].len() {
        let value: &str = match v["Thickness"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(FemError::ValueError),
        };
        let predicate: &str = match v["Thickness"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(FemError::PredicateError),
        };
        param.add_thickness_str(value, predicate);
    }

    for i in 0..v["YoungModulus"].len() {
        let value: &str = match v["YoungModulus"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(FemError::ValueError),
        };
        let predicate: &str = match v["YoungModulus"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(FemError::PredicateError),
        };
        param.add_young_modulus_str(value, predicate);
    }
    if mesh.is_2d() || mesh.is_3d() {

        for i in 0..v["PoissonRatio"].len() {
            let value: &str = match v["PoissonRatio"][i]["Value"].as_str() {
                Some(v) => v,
                None => return Err(FemError::ValueError),
            };
            let predicate: &str = match v["PoissonRatio"][i]["Predicate"].as_str() {
                Some(v) => v,
                None => return Err(FemError::PredicateError),
            };
            param.add_poisons_ratio_str(value, predicate);
        }
    }
    for i in 0..v["BoundaryConditions"].len() {
        let direct = get_json_direct("BoundaryConditions", &v, i)?;
        let value: &str = match v["BoundaryConditions"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(FemError::ValueError),
        };
        let predicate: &str = match v["BoundaryConditions"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(FemError::PredicateError),
        };
        param.add_boundary_condition_str(value, predicate, direct);
    }

    for i in 0..v["VolumeLoad"].len() {
        let direct = get_json_direct("VolumeLoad", &v, i)?;
        let value: &str = match v["VolumeLoad"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(FemError::ValueError),
        };
        let predicate: &str = match v["VolumeLoad"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(FemError::PredicateError),
        };
        param.add_volume_load_str(value, predicate, direct);
    }

    for i in 0..v["ConcentratedLoad"].len() {
        let direct = get_json_direct("ConcentratedLoad", &v, i)?;
        let value: &str = match v["ConcentratedLoad"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(FemError::ValueError),
        };
        let predicate: &str = match v["ConcentratedLoad"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(FemError::PredicateError),
        };
        param.add_concentrated_load_str(value, predicate, direct);
    }

    for i in 0..v["SurfaceLoad"].len() {
        let direct = get_json_direct("SurfaceLoad", &v, i)?;
        let value: &str = match v["SurfaceLoad"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(FemError::ValueError),
        };
        let predicate: &str = match v["SurfaceLoad"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(FemError::PredicateError),
        };
        param.add_surface_load_str(value, predicate, direct);
    }

    for i in 0..v["PressureLoad"].len() {
        let value: &str = match v["PressureLoad"][i]["Value"].as_str() {
            Some(v) => v,
            None => return Err(FemError::ValueError),
        };
        let predicate: &str = match v["PressureLoad"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(FemError::PredicateError),
        };
        param.add_pressure_load_str(value, predicate);
    }

    for i in 0..v["Variables"].len() {
        let s_val: &str = match v["Variables"][i].as_str() {
            Some(v) => v,
            None => return Err(FemError::ValueError),
        };
        let split: Vec<&str> = s_val.split(' ').collect();
        let name = split[0];
        let val = match split[1].to_string().parse() {
            Ok(v) => v,
            Err(_) => return Err(FemError::InvalidNumber),
        };
        param.add_variable(name, val);
    }

    for i in 0..v["StressStrainCurve"].len() {
        let mut value = Vec::<[f64; 2]>::new();
        for j in 0..v["StressStrainCurve"][i]["Value"].len() {
            let s = match v["StressStrainCurve"][i]["Value"][j][0].as_f64() {
                Some(v) => v,
                None => return Err(FemError::ValueError),
            };
            let e = match v["StressStrainCurve"][i]["Value"][j][1].as_f64() {
                Some(v) => v,
                None => return Err(FemError::ValueError),
            };
            value.push([s, e]);
        }
        let predicate: &str = match v["StressStrainCurve"][i]["Predicate"].as_str() {
            Some(v) => v,
            None => return Err(FemError::PredicateError),
        };
        param.add_stress_strain_curve_str(value, predicate);
        is_plasticity = true;
    }

    let res_name = match v["Result"].as_str() {
        Some(v) => v,
        None => return Err(FemError::ResultError),
    };

    if is_plasticity {
        generate::<FEMPlasticity>(&mesh, &param, res_name)
    } else {
        generate::<FEM>(&mesh, &param, res_name)
    }
}

#[allow(dead_code)]
fn get_json_direct<'a>(name: &'a str, v: &json::JsonValue, i: usize) -> Result<Direct, FemError> {
    let mut direct: Option<Direct> = None;
    let direct_str = match v[name][i]["Direct"].as_str() {
        Some(v) => v,
        None => return Err(FemError::DirectError),    
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
        None => Err(FemError::IncorrectDirectError),
    }
}

