// Информация о конечно-элементной сетке
//use std::str;
use std::fs::OpenOptions;
use std::io::{BufReader, BufWriter, prelude::*};
use ndarray::prelude::*;
use super::fe::FEType;
use super::error::FemError;
use super::msg::Messenger;
use super::util;

pub struct Mesh {
    pub freedom: usize,
    pub num_vertex: usize,
    pub num_fe: usize,
    pub num_be: usize,
    pub fe_type: FEType,
    pub x: Array2<f64>,
    pub fe: Array2<usize>,
    pub be: Array2<usize>,
    pub mesh_map: Vec<Vec<usize>>,
}

impl Mesh {
    pub fn new(file_name: &str) -> Result<Self, FemError> {
        let freedom: usize;
        let num_vertex: usize;
        let num_fe: usize;
        let mut num_be: usize;
        let fe_type;
        let fe_size;
        let be_size;
        let dim;
        let file = OpenOptions::new()
                .read(true)
                .write(true)
                .append(true)
                .create_new(false)
                .open(file_name)?;
        let mut reader = BufReader::new(&file);
        // Обработка типа КЭ
        let mut val = util::get_line(&mut reader)?;
        match val.trim() {
            "3" | "fe2d3" => {
                fe_type = FEType::FE2D3;
                dim = 2;
                fe_size = 3;
                be_size = 2;
                freedom = 2;
            }
            "4" | "fe3d4" => {
                fe_type = FEType::FE3D4;
                dim = 3;
                fe_size = 4;
                be_size = 3;
                freedom = 3;
            }
            "8" | "fe3d8" => {
                fe_type = FEType::FE3D8;
                dim = 3;
                fe_size = 8;
                be_size = 4;
                freedom = 3;
            }
            "24" | "fe2d4" => {
                fe_type = FEType::FE2D4;
                dim = 2;
                fe_size = 4;
                be_size = 2;
                freedom = 2;
            }
            "34" | "fe1d2" => {
                fe_type = FEType::FE1D2;
                dim = 1;
                fe_size = 2;
                be_size = 1;
                freedom = 1;
            }
            "223" | "fe3d3s" => {
                fe_type = FEType::FE3D3S;
                dim = 3;
                fe_size = 3;
                be_size = 0;
                freedom = 6;
            }
            "224" | "fe3d4s" => {
                fe_type = FEType::FE3D4S;
                dim = 3;
                fe_size = 4;
                be_size = 0;
                freedom = 6;
            }
            _ => return Err(FemError::InvalidFEType),
        }
        // Считывание количества узлов
        val = util::get_line(&mut reader)?;
        num_vertex = val.trim().parse()?;
        let mut x: Array2<f64> = Array2::zeros((num_vertex, dim));
        // Считывание узлов
        for i in 0..num_vertex {
            val = util::get_line(&mut reader)?;
            let ls = val.trim().split_whitespace();
            for (j, it) in ls.enumerate() {
                x[[i, j]] = it.parse()?;
            }
        }
        // Считывание количества КЭ
        val = util::get_line(&mut reader)?;
        num_fe = val.trim().parse()?;
        let mut fe: Array2<usize> = Array2::zeros((num_fe, fe_size));
        // Считывание КЭ
        for i in 0..num_fe {
            val = util::get_line(&mut reader)?;
            let ls = val.trim().split_whitespace();
            for (j, it) in ls.enumerate() {
                fe[[i, j]] = it.parse()?;
            }
        }
        // Считывание количества ГЭ
        val = util::get_line(&mut reader)?;
        num_be = val.trim().parse()?;
        let mut be: Array2<usize> = Array2::zeros((num_be, be_size));
        // Считывание ГЭ
        for i in 0..num_be {
            val = util::get_line(&mut reader)?;
            let ls = val.trim().split_whitespace();
            for (j, it) in ls.enumerate() {
                be[[i, j]] = it.parse()?;
            }
        }
        if fe_type == FEType::FE3D3S || fe_type == FEType::FE3D4S {
            num_be = num_fe;
            be = fe.clone();
        }
        // Считывание связей
        let mut mesh_map: Vec<Vec<usize>>;
        val = util::get_line(&mut reader)?;
        if val.trim().len() == 0 {
            // Если информации о связях в файле нет, создаем и записываем ее туда 
            let mut writer = BufWriter::new(file);
            mesh_map = Mesh::create_mesh_map(num_vertex, &fe);
            for row in &mesh_map {
                for elem in row {
                    writer.write(format!("{} ", elem).as_bytes())?;
                }
                writer.write(format!("\n").as_bytes())?;
            }
            if writer.flush().is_err() {
                return Err(FemError::WriteFile);
            }
        } else {
            mesh_map = vec![vec![]; num_vertex];
            for i in 0..num_vertex {
                if i > 0 {
                    val = util::get_line(&mut reader)?;
                }
                let ls = val.trim().split_whitespace();
                for (_, it) in ls.enumerate() {
                    mesh_map[i].push(it.parse()?);
                }
            }
        }

        println!("Mesh: {} ({}), nodes - {}, finite elements - {}", file_name, fe_type, num_vertex, num_fe);
        Ok(Self {freedom, num_vertex, num_fe, num_be, fe_type, x, fe, be, mesh_map})
    }
    // fn create_mesh_map(num_vertex: usize, fe: &Array2<usize>) -> Vec<Vec<usize>> {
    //     let mut link: Vec<bool> = vec![false; num_vertex];
    //     let mut mesh_map: Vec<Vec<usize>> = vec![vec![]; num_vertex];
    //     let mut msg = Messenger::new("Creating mesh map", 1, num_vertex as i64, 5);
    //     for i in 0..num_vertex {
    //         msg.add_progress();
    //         for j in 0..fe.shape()[0] {
    //             for k in 0..fe.shape()[1] {
    //                 if fe[[j, k]] == i {
    //                     for m in 0..fe.shape()[1] {
    //                         link[fe[[j, m]]] = true;
    //                     }
    //                 }
    //             }
    //         }
    //         for j in 0..link.len() {
    //             // if link[j] == true && j >= i {
    //             if link[j] == true {
    //                 mesh_map[i].push(j);
    //             }
    //         }
    //         link.fill(false);
    //     }
    //     // for i in 0..mesh_map.len() {
    //     //     for j in 0..mesh_map[i].len() {
    //     //         print!("{} ", mesh_map[i][j]);
    //     //     }            
    //     //     println!();
    //     // }
        
    //     // println!("{:?}", mesh_map);
    //     mesh_map
    // }
    fn create_mesh_map(num_vertex: usize, fe: &Array2<usize>) -> Vec<Vec<usize>> {
        let mut mesh_map: Vec<Vec<usize>> = vec![Vec::new(); num_vertex];
        let mut msg = Messenger::new("Creating mesh map", 1, fe.shape()[0] as i64, 5);
        for i in 0..fe.shape()[0] {
            msg.add_progress();
            for j in 0..fe.shape()[1] {
                for k in 0..fe.shape()[1] {
                    if mesh_map[fe[[i, j]]].iter().position(|l| *l == fe[[i, k]]) == None {
                        mesh_map[fe[[i, j]]].push(fe[[i, k]]);
                    }
                }
            }
        }
        for i in 0..num_vertex {
            mesh_map[i].sort();
        }
        //msg.stop();
        mesh_map
    }
    pub fn get_x_coord(&self, i: usize) -> Array1<f64> {
        let mut coord = Array1::zeros(self.x.shape()[1]);
        for j in 0..self.x.shape()[1] {
            coord[j] = self.x[[i, j]];
        }
        coord
    }
    pub fn get_fe_coord(&self, i: usize) -> Array2<f64> {
        let size1: usize = self.fe.shape()[1];
        let size2: usize = self.x.shape()[1];
        let mut x: Array2<f64> = Array2::zeros((size1, size2));
        for j in 0..size1 {
            for k in 0..size2 {
                x[[j, k]] = self.x[[self.fe[[i, j]] as usize, k]];
            }
        }
        x
    }
    pub fn get_be_coord(&self, i: usize) -> Array2<f64> {
        let size1: usize = self.be.shape()[1];
        let size2: usize = self.x.shape()[1];
        let mut x: Array2<f64> = Array2::zeros((size1, size2));
        for j in 0..size1 {
            for k in 0..size2 {
                x[[j, k]] = self.x[[self.be[[i, j]] as usize, k]];
            }
        }
        x
    }
    pub fn get_fe_center(&self, i: usize) -> Array1<f64> {
        let mut c;
        let mut coord = Array1::zeros(self.x.shape()[1]);
        for j in 0..self.x.shape()[1] {
            c = 0.;
            for k in 0..self.fe.shape()[1] {
                c += self.x[[self.fe[[i, k]], j]];
            }
            coord[j] = c / self.fe.shape()[1] as f64;
        }
        coord
    }
    pub fn get_be_center(&self, i: usize) -> Array1<f64> {
        let mut c;
        let mut coord = Array1::zeros(self.x.shape()[1]);
        for j in 0..self.x.shape()[1] {
            c = 0.;
            for k in 0..self.be.shape()[1] {
                c += self.x[[self.be[[i, k]], j]];
            }
            coord[j] = c / self.be.shape()[1] as f64;
        }
        coord
    }
    pub fn is_1d(&self) -> bool {
        if self.fe_type == FEType::FE1D2 { true } else { false }
    }
    pub fn is_2d(&self) -> bool {
        if self.fe_type == FEType::FE2D3 || self.fe_type == FEType::FE2D4 { true } else { false }
    }
    pub fn is_3d(&self) -> bool {
        if self.fe_type == FEType::FE3D3S || self.fe_type == FEType::FE3D4S || self.fe_type == FEType::FE3D4 || self.fe_type == FEType::FE3D8 { true } else { false }
    }
    pub fn is_shell(&self) -> bool {
        if self.fe_type == FEType::FE3D3S || self.fe_type == FEType::FE3D4S { true } else { false }
    }
    pub fn fe_volume(&self, i: usize) -> f64 {
        let x = self.get_fe_coord(i);
        match self.fe_type {
            FEType::FE1D2 => self.volume_1d2(x),
            FEType::FE2D3 | FEType::FE3D3S => self.volume_2d3(x),
            FEType::FE2D4 | FEType::FE3D4S => self.volume_2d4(x),
            FEType::FE3D4 => self.volume_3d4(x),
            FEType::FE3D8 => self.volume_3d8(x),
        }
    }
    // Длина отрезка
    pub fn volume_1d2(&self, x: Array2<f64>) -> f64 {
        ((x[[0, 0]] - x[[1, 0]]) * (x[[0, 0]] - x[[1, 0]])).sqrt()
    }
    // Площадь треугольника
    pub fn volume_2d3(&self, x: Array2<f64>) -> f64 {
        let a = ((x[[0, 0]] - x[[1, 0]]) * (x[[0, 0]] - x[[1, 0]]) + (x[[0, 1]] - x[[1, 1]]) * (x[[0, 1]] - x[[1, 1]]) + if x.shape()[1] == 3 { (x[[0, 2]] - x[[1, 2]]) * (x[[0, 2]] - x[[1, 2]]) } else { 0. }).sqrt();
        let b = ((x[[0, 0]] - x[[2, 0]]) * (x[[0, 0]] - x[[2, 0]]) + (x[[0, 1]] - x[[2, 1]]) * (x[[0, 1]] - x[[2, 1]]) + if x.shape()[1] == 3 { (x[[0, 2]] - x[[2, 2]]) * (x[[0, 2]] - x[[2, 2]]) } else { 0. }).sqrt();
        let c = ((x[[2, 0]] - x[[1, 0]]) * (x[[2, 0]] - x[[1, 0]]) + (x[[2, 1]] - x[[1, 1]]) * (x[[2, 1]] - x[[1, 1]]) + if x.shape()[1] == 3 { (x[[2, 2]] - x[[1, 2]]) * (x[[2, 2]] - x[[1, 2]]) } else { 0. }).sqrt();
        // let a = ((x[[0, 0]] - x[[1, 0]]) * (x[[0, 0]] - x[[1, 0]]) + (x[[0, 1]] - x[[1, 1]]) * (x[[0, 1]] - x[[1, 1]])).sqrt();
        // let b = ((x[[0, 0]] - x[[2, 0]]) * (x[[0, 0]] - x[[2, 0]]) + (x[[0, 1]] - x[[2, 1]]) * (x[[0, 1]] - x[[2, 1]])).sqrt();
        // let c = ((x[[2, 0]] - x[[1, 0]]) * (x[[2, 0]] - x[[1, 0]]) + (x[[2, 1]] - x[[1, 1]]) * (x[[2, 1]] - x[[1, 1]])).sqrt();
        let p = 0.5 * (a + b + c);
        (p * (p - a) * (p - b) * (p - c)).sqrt()
    }
    // Площадь четырехугольника
    pub fn volume_2d4(&self, x: Array2<f64>) -> f64 {
        self.volume_2d3(array![ [x[[0, 0]], x[[0, 1]], if x.shape()[1] == 3 { x[[0, 2]] } else { 0. } ], [x[[1, 0]], x[[1, 1]], if x.shape()[1] == 3 { x[[1, 2]] } else { 0. }], [x[[2, 0]], x[[2, 1]], if x.shape()[1] == 3 { x[[2, 2]] } else { 0. } ] ]) + 
        self.volume_2d3(array![ [x[[2, 0]], x[[2, 1]], if x.shape()[1] == 3 { x[[2, 2]] } else { 0. } ], [x[[3, 0]], x[[3, 1]], if x.shape()[1] == 3 { x[[3, 2]] } else { 0. }], [x[[0, 0]], x[[0, 1]], if x.shape()[1] == 3 { x[[0, 2]] } else { 0. } ] ])
    }
    // Объем пирамиды
    pub fn volume_3d4(&self, x: Array2<f64>) -> f64 {
        let m: Array2<f64> = array![ [ x[[1, 0]] - x[[0, 0]], x[[1, 1]] - x[[0, 1]], x[[1, 2]] - x[[0, 2]] ],
        [ x[[2, 0]] - x[[0, 0]], x[[2, 1]] - x[[0, 1]], x[[2, 2]] - x[[0, 2]] ],
        [ x[[3, 0]] - x[[0, 0]], x[[3, 1]] - x[[0, 1]], x[[3, 2]] - x[[0, 2]] ] ];
        util::det(&m).unwrap().abs() / 6.0
    }
    // Объем четырехугольного шестигранника
    pub fn volume_3d8(&self, x: Array2<f64>) -> f64 {
        let ind: Array2<usize> = array![ [0, 1, 4, 7], [4, 1, 5, 7], [1, 2, 6, 7], [1, 5, 6, 7], [1, 2, 3, 7], [0, 3, 1, 7] ];
        let mut v = 0.;
        for i in 0..6 {
            let m: Array2<f64> = array![ [x[[ind[[i, 0]], 0]], x[[ind[[i, 0]], 1]], x[[ind[[i, 0]], 2]]], 
                [x[[ind[[i, 1]], 0]], x[[ind[[i, 1]], 1]], x[[ind[[i, 1]], 2]]], 
                [x[[ind[[i, 2]], 0]], x[[ind[[i, 2]], 1]], x[[ind[[i, 2]], 2]]], 
                [x[[ind[[i, 3]], 0]], x[[ind[[i, 3]], 1]], x[[ind[[i, 3]], 2]]] ];
            v += self.volume_3d4(m);
        }
        v
    }
    pub fn be_volume(&self, i: usize) -> f64 {
        let x = self.get_be_coord(i);
        match self.fe_type {
            FEType::FE1D2 => 1.0,
            FEType::FE2D3 | FEType::FE2D4 => self.volume_1d2(x),
            FEType::FE3D3S | FEType::FE3D4 => self.volume_2d3(x),
            FEType::FE3D4S | FEType::FE3D8 => self.volume_2d4(x),
        }
    }
    pub fn be_normal(&self, i: usize) -> Array1<f64>{
        let x = self.get_be_coord(i);
        let normal: Array1<f64> = match self.fe_type {
            FEType::FE1D2 => array![1., 0., 0.],
            FEType::FE2D3 | FEType::FE2D4 => 
                array![self.x[[0, 1]] - self.x[[1, 1]], self.x[[1, 0]] - self.x[[0, 0]], 0.],
            FEType::FE3D3S | FEType::FE3D4S | FEType::FE3D4 | FEType::FE3D8 => 
                array![ (x[[1, 1]] - x[[0, 1]]) * (x[[2, 2]] - x[[0, 2]]) - (x[[2, 1]] - x[[0, 1]]) * (x[[1, 2]] - x[[0, 2]]), 
                    (x[[2, 0]] - x[[0, 0]]) * (x[[1, 2]] - x[[0, 2]]) - (x[[1, 0]] - x[[0, 0]]) * (x[[2, 2]] - x[[0, 2]]), 
                    (x[[1, 0]] - x[[0, 0]]) * (x[[2, 1]] - x[[0, 1]]) - (x[[2, 0]] - x[[0, 0]]) * (x[[1, 1]] - x[[0, 1]]) ],
        };
        let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();   
        normal / len  
    }
    pub fn nnz(&self) -> usize {
        let mut val = 0;
        for i in &self.mesh_map {
            val += i.len();
        }
        val * self.freedom * self.freedom
    }
}


// fn get_line<R: BufRead>(reader: &mut R) -> Result<String, Error> {
//     let mut res;
//     let mut bytes: Vec<u8> = Vec::new();
//     loop {
//         let len = match reader.read_until(b'\n', &mut bytes) {
//             Err(_) => return Err(Error::ReadFile),
//             Ok(len) => len,
//         };
//         if len == 0 {
//             return Ok(String::new());
//         }
//         res = match str::from_utf8(&bytes) {
//             Ok(res) => res,
//             Err(_) => return Err(Error::InvalidNumber),
//         };
//         if res.trim().len() > 0 {
//             break;
//         }
//     }
//     Ok(res.to_string())
// }
