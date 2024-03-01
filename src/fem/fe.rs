// Поддержка конечного элемента
use std::fmt;
use ndarray::prelude::*;
use super::util;
use super::error::FemError;

// Типы конечных элементов (КЭ)
#[derive(Copy, Clone, PartialEq)]
pub enum FEType {
    FE1D2,
    FE2D3,
    FE2D4,
    FE3D3S,
    FE3D4S,
    FE3D4,
    FE3D8,
}

impl fmt::Display for FEType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let str = match self {
            FEType::FE1D2 => "fe1d2",     
            FEType::FE2D3 => "fe2d3",     
            FEType::FE2D4 => "fe2d4",     
            FEType::FE3D3S => "fe3d3s",     
            FEType::FE3D4S => "fe3d4s",     
            FEType::FE3D4 => "fe3d4",     
            FEType::FE3D8 => "fe3d8",     
        };
        write!(f, "{}", str)
    }
}

pub trait FiniteElement {
    fn size(&self) -> usize;
    fn freedom(&self) -> usize;
    fn shape_factor(&self, i: usize, j: usize) -> f64;
    fn x(&self, i: usize, j: usize) -> f64;
    fn w(&self) -> Array1<f64>;
    fn e(&self) -> Array1<f64>;
    fn create(&self) -> Result<Array2<f64>, FemError> {
        let mut a: Array2<f64> = Array2::zeros((self.size(), self.size()));
        let mut b: Array1<f64> = Array1::zeros(self.size());
        let mut c: Array2<f64> = Array2::zeros((self.size(), self.size()));
        for i in 0..self.size() {
            for j in 0..self.size() {
                for k in 0..self.size() {
                    a[[j, k]] = self.shape_factor(j, k);
                }
                b[j] = if i == j { 1. } else { 0. };
            }
            let x = util::solve(&mut a, &mut b, 1.0e-10)?;
            for j in 0..self.size() {
                c[[j, i]] = x[j];
            }
        }
        Ok(c)
    }
} 


pub mod fe1d {
    use super::FemError; 
    use super::FiniteElement;
    use ndarray::prelude::*;
    use crate::fem::util;

    pub struct FE1D2 {
        e: [f64; 2],        // Модуль Юнга и коэффициент Пуассона
        thk: f64,           // Толщина
        x: Array2<f64>,     // Координаты вершин
    }
    
    pub trait FiniteElement1D: FiniteElement {
        fn thk(&self) -> f64;
        fn d(&self) -> Array2<f64> {
            array![[self.e()[0]]]      
        }
        fn jacobi(&self, _i: usize) -> Array2<f64> {
            array![[(self.x(1, 0) - self.x(0, 0)) * 0.5]]
        }
        fn dx(&self, c: &Array2<f64>, i: usize) -> f64;
        fn dxi(&self, i: usize, j: usize) -> f64;
        fn generate(&mut self) -> Result<Array2<f64>, FemError> {
            // self.create()?;
            let size = self.size() * self.freedom();
            let mut local: Array2<f64> = Array2::zeros((size , size));
            for i in 0..self.w().len() {
                let mut b: Array2<f64> = Array2::zeros((1, size));
                let jacobi = self.jacobi(i);
                let inv_jacobi = util::inv(&jacobi)?;
                let jacobian = util::det(&jacobi)?;
                for j in 0..self.size() {
                    b[[0, j]] = inv_jacobi[[0, 0]] * self.dxi(i, j);
                }
                local = local + b.t().dot(&self.d()).dot(&b) * (self.w()[i] * self.thk() * jacobian.abs());
            }
            Ok(local)
        }
        fn calc(&self, u: &Array1<f64>) -> Result<Array2<f64>, FemError> {
            let c = self.create()?;
            // println!("\n{:?}", c);
            let mut res: Array2<f64> = Array2::zeros((2, self.size()));
            for i in 0..self.size() {
                let mut b: Array2<f64> = Array2::zeros((1, self.size())); 
                for j in 0..self.size() {
                    b[[0, j]] = self.dx(&c, j);
                }
                // println!("{:?}", b);
                let e = b.dot(u);
                let s = self.d().dot(&e);
                // println!("{:?}", e);
                // println!("{:?}", s);
                res[[0, i]] += e[0];
                res[[1, i]] += s[0];
            }
            Ok(res)
        }
    }
    
    impl FiniteElement for FE1D2 {
        fn w(&self) -> Array1<f64> {
            array![ 0.55555555556, 0.88888888889, 0.55555555556 ]    
        }
        fn e(&self) -> Array1<f64> {
            array![ self.e[0], self.e[1] ]
        }
        fn x(&self, i: usize, _j: usize) -> f64 {
            self.x[[i, 0]]
        }
        fn size(&self) -> usize {
            2
        }
        fn freedom(&self) -> usize {
            1
        }
        fn shape_factor(&self, i: usize, j: usize) -> f64 {
            [ 1.0, self.x[[i, 0]] ][j]
        }
    }
    
    impl FiniteElement1D for FE1D2 {
        fn thk(&self) -> f64 {
            self.thk
        }
        fn dxi(&self, _i: usize, j: usize) -> f64 {
            [ -0.5, 0.5 ][j]
        }
        fn dx(&self, c: &Array2<f64>, j: usize) -> f64 {
            c[[1, j]]
        }
    
    }
    impl FE1D2 {
        pub fn new(x: Array2<f64>, e: [f64; 2], thk: f64) -> Self {
            Self { x, e, thk, }
        }
    }
}

pub mod fe2d {
    use super::FemError; 
    use super::FiniteElement;
    use ndarray::prelude::*;
    use crate::fem::util;

    pub struct FE2D3 {
        is_shell: bool,
        thk: f64,           
        e: [f64; 2],        
        x: Array2<f64>,   
        m: Array2<f64>,     
    }
    
    pub struct FE2D4 {
        is_shell: bool,
        thk: f64,           
        e: [f64; 2],        
        x: Array2<f64>,   
        m: Array2<f64>,     
    }
    
    pub trait FiniteElement2D: FiniteElement {
        fn is_shell(&self) -> bool;
        fn thk(&self) -> f64;
        fn shape(&self, i: usize, j: usize) -> f64;
        fn m(&self) -> &Array2<f64>;
        fn ed(&self) -> Array2<f64> {
            array![[ self.e()[0] / (2. + 2. * self.e()[1]),  0. ],
                   [ 0.,  self.e()[0] / (2. + 2. * self.e()[1]) ]]
        }
        fn dx(&self, c: &Array2<f64>, i: usize, j: usize) -> f64;
        fn dy(&self, c: &Array2<f64>, i: usize, j: usize) -> f64;
        fn dxi(&self, i: usize, j: usize) -> f64;
        fn deta(&self, i: usize, j: usize) -> f64;
        fn d(&self) -> Array2<f64> {
            array![[ self.e()[0] / (1.0 - self.e()[1] * self.e()[1]),  self.e()[1] * self.e()[0] / (1.0 - self.e()[1] * self.e()[1]), 0. ],
                [ self.e()[1] * self.e()[0] / (1.0 - self.e()[1] * self.e()[1]), self.e()[0] / (1.0 - self.e()[1] * self.e()[1]),  0. ], 
                [ 0., 0., 0.5 * (1.0 - self.e()[1]) * self.e()[0] / (1.0 - self.e()[1] * self.e()[1])]]
        }
        fn jacobi(&self, i: usize) -> Array2<f64> {
            let mut res: Array2<f64> = Array2::zeros((2, 2));
            for j in 0..2 {
                for k in 0..self.size() {
                    res[[0, j]] += self.dxi(i, k) * self.x(k, j);
                    res[[1, j]] += self.deta(i, k) * self.x(k, j);
                }
            }
            res
        }
        fn generate(&mut self) -> Result<Array2<f64>, FemError> {
            let size = self.size() * self.freedom();
            let mut local: Array2<f64> = Array2::zeros((size , size));
            for i in 0..self.w().len() {
                // Якобиан и обратная матрица Якоби
                let jacobi = self.jacobi(i);
                let inv_jacobi = util::inv(&jacobi)?;
                let jacobian = util::det(&jacobi)?;
                if !self.is_shell() {
                    let mut b: Array2<f64> = Array2::zeros((3, size));
                    for j in 0..self.size() {
                        b[[0, j * self.freedom() + 0]] = inv_jacobi[[0, 0]] * self.dxi(i, j) + inv_jacobi[[0, 1]] * self.deta(i, j);
                        b[[1, j * self.freedom() + 1]] = inv_jacobi[[1, 0]] * self.dxi(i, j) + inv_jacobi[[1, 1]] * self.deta(i, j);
                        b[[2, j * self.freedom() + 1]] = b[[0, j * self.freedom() + 0]];
                        b[[2, j * self.freedom() + 0]] = b[[1, j * self.freedom() + 1]];
                    }
                    local = local + b.t().dot(&self.d()).dot(&b) * (self.w()[i] * self.thk() * jacobian.abs());
                } else {
                    let mut bm: Array2<f64> = Array2::zeros((3, size));
                    let mut bp: Array2<f64> = Array2::zeros((3, size));
                    let mut bc: Array2<f64> = Array2::zeros((2, size));
                    for j in 0..self.size()  {
                        bm[[0, self.freedom() * j + 0]] = inv_jacobi[[0, 0]] * self.dxi(i, j) + inv_jacobi[[0, 1]] * self.deta(i, j);
                        bm[[2, self.freedom() * j + 1]] = inv_jacobi[[0, 0]] * self.dxi(i, j) + inv_jacobi[[0, 1]] * self.deta(i, j);
                        bp[[0, self.freedom() * j + 3]] = inv_jacobi[[0, 0]] * self.dxi(i, j) + inv_jacobi[[0, 1]] * self.deta(i, j);
                        bp[[2, self.freedom() * j + 4]] = inv_jacobi[[0, 0]] * self.dxi(i, j) + inv_jacobi[[0, 1]] * self.deta(i, j);
                        bc[[0, self.freedom() * j + 2]] = inv_jacobi[[0, 0]] * self.dxi(i, j) + inv_jacobi[[0, 1]] * self.deta(i, j);
                        bm[[1, self.freedom() * j + 1]] = inv_jacobi[[1, 0]] * self.dxi(i, j) + inv_jacobi[[1, 1]] * self.deta(i, j);
                        bm[[2, self.freedom() * j + 0]] = inv_jacobi[[1, 0]] * self.dxi(i, j) + inv_jacobi[[1, 1]] * self.deta(i, j);
                        bp[[1, self.freedom() * j + 4]] = inv_jacobi[[1, 0]] * self.dxi(i, j) + inv_jacobi[[1, 1]] * self.deta(i, j);
                        bp[[2, self.freedom() * j + 3]] = inv_jacobi[[1, 0]] * self.dxi(i, j) + inv_jacobi[[1, 1]] * self.deta(i, j);
                        bc[[1, self.freedom() * j + 2]] = inv_jacobi[[1, 0]] * self.dxi(i, j) + inv_jacobi[[1, 1]] * self.deta(i, j);
                        bc[[0, self.freedom() * j + 3]] = self.shape(i, j);
                        bc[[1, self.freedom() * j + 4]] = self.shape(i, j);
                    }
                    local = local + ((bm.t().dot(&self.d()).dot(&bm)) * self.thk() + 
                                    (bp.t().dot(&self.d()).dot(&bp)) * (self.thk().powf(3.) / 12.) + 
                                    (bc.t().dot(&self.ed()).dot(&bc)) * (self.thk() * 5. / 6.)) * (self.w()[i] * jacobian.abs());
                }
            }
            if self.is_shell() {
                // Поиск максимального диагонального элемента
                let mut singular = 0.;
                for i in 0..size {
                    if local[[i, i]] > singular {
                        singular = local[[i, i]];
                    }
                }
                singular *= 1.0E-3;
                // Устранение сингулярности
                for i in 0..self.size() {
                    local[[self.freedom() * (i + 1) - 1, self.freedom() * (i + 1) - 1]] = singular;
                }
                // Преобразование из локальных координат в глобальные
                let m = util::create_ext_transform_matrix(&self.m(), self.size(), self.freedom());
                local = m.t().dot(&local).dot(&m);
            }
            // println!("\n{:?}", local);
            Ok(local)
        }
        fn calc(&self, u: &Array1<f64>) -> Result<Array2<f64>, FemError> {
            let c = self.create()?;
            let mut res: Array2<f64> = Array2::zeros((if self.is_shell() { 12 } else { 6 }, self.size() * self.freedom()));
            if !self.is_shell() {
                for i in 0..self.size() {
                    let mut b: Array2<f64> = Array2::zeros((3, self.size() * self.freedom())); 
                    for j in 0..self.size() {
                        b[[0, j * self.freedom() + 0]] = self.dx(&c, i, j);
                        b[[1, j * self.freedom() + 1]] = self.dy(&c, i, j);
                        b[[2, j * self.freedom() + 0]] = self.dy(&c, i, j);
                        b[[2, j * self.freedom() + 1]] = self.dx(&c, i, j);
                    }
                    let e = b.dot(u);
                    let s = self.d().dot(&e);
                    for j in 0..3 {
                        res[[j, i]] += e[j];
                        res[[j + 3, i]] += s[j];
                    }
                }
            } else {
                // Преобразование перемещений
                let m = util::create_ext_transform_matrix(&self.m(), self.size(), self.freedom());
                let lu = m.dot(u);    
                for i in 0..self.size() {
                    let mut bm = Array2::<f64>::zeros((3, self.size() * self.freedom()));
                    let mut bp = Array2::<f64>::zeros((3, self.size() * self.freedom()));
                    let mut bc = Array2::<f64>::zeros((2, self.size() * self.freedom()));
                    for j in 0..self.size() {
                        bm[[0, self.freedom() * j + 0]] = self.dx(&c, i, j);
                        bm[[2, self.freedom() * j + 1]] = self.dx(&c, i, j);
                        bp[[2, self.freedom() * j + 4]] = self.dx(&c, i, j);
                        bp[[0, self.freedom() * j + 3]] = self.dx(&c, i, j);
                        bc[[0, self.freedom() * j + 2]] = self.dx(&c, i, j);
                        bm[[1, self.freedom() * j + 1]] = self.dy(&c, i, j);
                        bm[[2, self.freedom() * j + 0]] = self.dy(&c, i, j);
                        bp[[1, self.freedom() * j + 4]] = self.dy(&c, i, j);
                        bp[[2, self.freedom() * j + 3]] = self.dy(&c, i, j);
                        bc[[1, self.freedom() * j + 2]] = self.dy(&c, i, j);
                        bc[[0, self.freedom() * j + 3]] = if i == j { 1. } else { 0. };
                        bc[[1, self.freedom() * j + 4]] = if i == j { 1. } else { 0. };
                    }
                    let strainm = bm.dot(&lu);
                    let strainp = bp.dot(&lu);
                    let strainc = bc.dot(&lu);
                    let stressm = self.d().dot(&strainm);
                    let stressp = self.d().dot(&strainp) * self.thk() * 0.5;
                    let stressc = self.ed().dot(&strainc);
                    let local_strain = array![[strainm[0] + strainp[0], strainm[2] + strainp[2], strainc[0]],
                                              [strainm[2] + strainp[2], strainm[1] + strainp[1], strainc[1]],
                                              [strainc[0], strainc[1], 0.]];
                    let local_stress = array![[stressm[0] + stressp[0], stressm[2] + stressp[2], stressc[0]],
                                              [stressm[2] + stressp[2], stressm[1] + stressp[1], stressc[1]],
                                              [stressc[0], stressc[1], 0.]];
                    let global_strain = self.m().t().dot(&local_strain).dot(self.m());
                    let global_stress = self.m().t().dot(&local_stress).dot(self.m());
                    let index = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
                    for j in 0..6 {
                        res[[j, i]] += global_strain[[index[j].0, index[j].1]]; 
                        res[[6 + j, i]] += global_stress[[index[j].0, index[j].1]]; 
                    }

                    // res[[0, i]] += global_strain[[0, 0]];    // Exx
                    // res[[1, i]] += global_strain[[1, 1]];    // Eyy
                    // res[[2, i]] += global_strain[[2, 2]];    // Ezz
                    // res[[3, i]] += global_strain[[0, 1]];    // Exy
                    // res[[4, i]] += global_strain[[0, 2]];    // Exz
                    // res[[5, i]] += global_strain[[1, 2]];    // Eyz
                    // res[[6, i]] += global_stress[[0, 0]];    // Sxx
                    // res[[7, i]] += global_stress[[1, 1]];    // Syy
                    // res[[8, i]] += global_stress[[2, 2]];    // Szz
                    // res[[9, i]] += global_stress[[0, 1]];    // Sxy
                    // res[[10, i]] += global_stress[[0, 2]];   // Sxz
                    // res[[11, i]] += global_stress[[1, 2]];   // Syz
                }
            }

            Ok(res)
        }
    }

    impl FiniteElement for FE2D4 {
        fn w(&self) -> Array1<f64> {
            array![ 1.0, 1.0, 1.0, 1.0  ]  
        }
        fn e(&self) -> Array1<f64> {
            array![ self.e[0], self.e[1] ]
        }
        fn x(&self, i: usize, j: usize) -> f64 {
            self.x[[i, j]]
        }
        fn size(&self) -> usize {
            4
        }
        fn freedom(&self) -> usize {
            if !self.is_shell { 2 } else { 6 }
        }
        fn shape_factor(&self, i: usize, j: usize) -> f64 {
            [ 1.0, self.x[[i, 0]], self.x[[i, 1]], self.x[[i, 0] ] * self.x[[i, 1]]][j]
        }
    }

    impl FiniteElement for FE2D3 {
        fn w(&self) -> Array1<f64> {
            array![ 0.166666666667, 0.166666666667, 0.166666666667 ]  
        }
        fn e(&self) -> Array1<f64> {
            array![ self.e[0], self.e[1] ]
        }
        fn x(&self, i: usize, j: usize) -> f64 {
            self.x[[i, j]]
        }
        fn size(&self) -> usize {
            3
        }
        fn freedom(&self) -> usize {
            if !self.is_shell { 2 } else { 6 }
        }
        fn shape_factor(&self, i: usize, j: usize) -> f64 {
            [ 1.0, self.x[[i, 0]], self.x[[i, 1]] ][j]
        }
    }
    
    impl FiniteElement2D for FE2D3 {
        fn is_shell(&self) -> bool {
            self.is_shell
        }
        fn thk(&self) -> f64 {
            self.thk
        }
        fn m(&self) -> &Array2<f64> {
            &self.m
        }
        fn dx(&self, c: &Array2<f64>, _i: usize, j: usize) -> f64 {
            c[[1, j]]
        }
        fn dy(&self, c: &Array2<f64>, _i: usize, j: usize) -> f64 {
            c[[2, j]]
        }
        fn dxi(&self, _i: usize, j: usize) -> f64 {
            [ -1.0, 1.0, 0.0 ][j]
        }
        fn deta(&self, _i: usize, j: usize) -> f64 {
            [ -1.0, 0.0, 1.0 ][j]
        }
        fn shape(&self, i: usize, j: usize) -> f64 {
            let xi = array![ 0.0, 0.5, 0.5 ];
            let eta = array![ 0.5, 0.0, 0.5 ];
            [ 1.0 - xi[i] - eta[i], xi[i], eta[i] ][j]
        }
    }
    
    impl FiniteElement2D for FE2D4 {
        fn is_shell(&self) -> bool {
            self.is_shell
        }
        fn thk(&self) -> f64 {
            self.thk
        }
        fn m(&self) -> &Array2<f64> {
            &self.m
        }
        fn dx(&self, c: &Array2<f64>, i: usize, j: usize) -> f64 {
            c[[1, j]] + c[[3, j]] * self.x[[i, 1]]
        }
        fn dy(&self, c: &Array2<f64>, i: usize, j: usize) -> f64 {
            c[[2, j]] + c[[3, j]] * self.x[[i, 0]]
        }
        fn dxi(&self, i: usize, j: usize) -> f64 {
            let eta = array![ -0.57735027, 0.57735027, -0.57735027, 0.57735027 ];
            [ -0.25 * (1.0 - eta[i]), 0.25 * (1.0 - eta[i]), 0.25 * (1.0 + eta[i]), -0.25 * (1.0 + eta[i]) ][j]
        }
        fn deta(&self, i: usize, j: usize) -> f64 {
            let xi = array![ -0.57735027, -0.57735027, 0.57735027, 0.57735027 ];
            [ -0.25 * (1.0 - xi[i]), -0.25 * (1.0 + xi[i]), 0.25 * (1.0 + xi[i]), 0.25 * (1.0 - xi[i]) ][j]
        }
        fn shape(&self, i: usize, j: usize) -> f64 {
            let xi = array![ -0.57735027, -0.57735027, 0.57735027, 0.57735027 ];
            let eta = array![ -0.57735027, 0.57735027, -0.57735027, 0.57735027 ];
            [ 0.25 * (1.0 - xi[i]) * (1.0 - eta[i]), 0.25 * (1.0 + xi[i]) * (1.0 - eta[i]), 0.25 * (1.0 + xi[i]) * (1.0 + eta[i]), 0.25 * (1.0 - xi[i]) * (1.0 + eta[i]) ][j]
        }
    }

    impl FE2D3 {
        pub fn new(x: Array2<f64>, e: [f64; 2], thk: f64, is_shell: bool, ) -> Self {
            let mut m = Array2::<f64>::zeros((0, 0));
            let mut x = x;
            if is_shell {
                m = util::create_transform_matrix(&x);
                let y = m.dot(&x.t());
                for i in 0..x.shape()[0] {
                    for j in 0..x.shape()[1] {
                        x[[i, j]] = y[[j, i]];
                    }
                }
            }
            Self { x, e, thk, is_shell, m }
        }
    }
    
    impl FE2D4 {
        pub fn new(x: Array2<f64>, e: [f64; 2], thk: f64, is_shell: bool, ) -> Self {
            let mut m = Array2::<f64>::zeros((0, 0));
            let mut x = x;
            if is_shell {
                m = util::create_transform_matrix(&x);
                let y = m.dot(&x.t());
                for i in 0..x.shape()[0] {
                    for j in 0..x.shape()[1] {
                        x[[i, j]] = y[[j, i]];
                    }
                }
            }
            Self { x, e, thk, is_shell, m }
        }
    }
}

pub mod fe3d {
    use super::FemError; 
    use super::FiniteElement;
    use ndarray::prelude::*;
    use crate::fem::util;

    pub struct FE3D4 {
        e: [f64; 2],        
        x: Array2<f64>,     
    }
    
    pub struct FE3D8 {
        e: [f64; 2],        
        x: Array2<f64>,     
    }

    pub trait FiniteElement3D: FiniteElement {
        fn dx(&self, c: &Array2<f64>, i: usize, j: usize) -> f64;
        fn dy(&self, c: &Array2<f64>, i: usize, j: usize) -> f64;
        fn dz(&self, c: &Array2<f64>, i: usize, j: usize) -> f64;
        fn dxi(&self, i: usize, j: usize) -> f64;
        fn deta(&self, i: usize, j: usize) -> f64;
        fn dpsi(&self, i: usize, j: usize) -> f64;
        fn jacobi(&self, i: usize) -> Array2<f64> {
            let mut res: Array2<f64> = Array2::zeros((3, 3));
            for j in 0..3 {
                for k in 0..self.size() {
                    res[[0, j]] += self.dxi(i, k) * self.x(k, j);
                    res[[1, j]] += self.deta(i, k) * self.x(k, j);
                    res[[2, j]] += self.dpsi(i, k) * self.x(k, j);
                }
            }
            res
        }
        fn d(&self) -> Array2<f64> {
            let e = self.e()[0];
            let m = self.e()[1];
            array![[ e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), m / (1.0 - m) * e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), m / (1.0 - m) * e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), 0.0, 0.0, 0.0 ], 
                [ m / (1.0 - m) * e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), m / (1.0 - m) * e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), 0.0, 0.0, 0.0 ],
                [ m / (1.0 - m) * e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), m / (1.0 - m) * e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), 0.0, 0.0, 0.0 ],
                [ 0.0, 0.0, 0.0, 0.5 * (1.0 - 2.0 * m) / (1.0 - m) * e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), 0.0, 0.0 ],
                [ 0.0, 0.0, 0.0, 0.0, 0.5 * (1.0 - 2.0 * m) / (1.0 - m) * e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m), 0.0 ],
                [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.5 * (1.0 - 2.0 * m) / (1.0 - m) * e * (1.0 - m) / (1.0 + m) / (1.0 - 2.0 * m)]]
        }
        fn generate(&mut self) -> Result<Array2<f64>, FemError> {
            // self.create()?;
            let size = self.size() * self.freedom();
            let mut local: Array2<f64> = Array2::zeros((size , size));
            for i in 0..self.w().len() {
                let mut b: Array2<f64> = Array2::zeros((6, size));
                // Якобиан и обратная матрица Якоби
                let jacobi = self.jacobi(i);
                let inv_jacobi = util::inv(&jacobi)?;
                let jacobian = util::det(&jacobi)?;
                for j in 0..self.size()  {
                    b[[0, j * self.freedom() + 0]] = inv_jacobi[[0, 0]] * self.dxi(i, j) + inv_jacobi[[0, 1]] * self.deta(i, j) + inv_jacobi[[0, 2]] * self.dpsi(i, j);
                    b[[3, j * self.freedom() + 1]] = b[[0, j * self.freedom() + 0]];
                    b[[5, j * self.freedom() + 2]] = b[[0, j * self.freedom() + 0]];
                    b[[1, j * self.freedom() + 1]] = inv_jacobi[[1, 0]] * self.dxi(i, j) + inv_jacobi[[1, 1]] * self.deta(i, j) + inv_jacobi[[1, 2]] * self.dpsi(i, j);
                    b[[3, j * self.freedom() + 0]] = b[[1, j * self.freedom() + 1]];
                    b[[4, j * self.freedom() + 2]] = b[[1, j * self.freedom() + 1]];
                    b[[2, j * self.freedom() + 2]] = inv_jacobi[[2, 0]] * self.dxi(i, j) + inv_jacobi[[2, 1]] * self.deta(i, j) + inv_jacobi[[2, 2]] * self.dpsi(i, j);
                    b[[4, j * self.freedom() + 1]] = b[[2, j * self.freedom() + 2]];
                    b[[5, j * self.freedom() + 0]] = b[[2, j * self.freedom() + 2]];
                }
                local = local + b.t().dot(&self.d()).dot(&b) * (self.w()[i] * jacobian.abs());
            }
            Ok(local)
        }
        fn calc(&self, u: &Array1<f64>) -> Result<Array2<f64>, FemError> {
            let c = self.create()?;
            let mut res: Array2<f64> = Array2::zeros((12, self.size() * self.freedom()));
            for i in 0..self.size() {
                let mut b: Array2<f64> = Array2::zeros((6, self.size() * self.freedom())); 
                for j in 0..self.size() {
                    b[[0, j * self.freedom() + 0]] = self.dx(&c, i, j);
                    b[[3, j * self.freedom() + 0]] = self.dy(&c, i, j);
                    b[[5, j * self.freedom() + 0]] = self.dz(&c, i, j);
                    b[[3, j * self.freedom() + 1]] = self.dx(&c, i, j);
                    b[[1, j * self.freedom() + 1]] = self.dy(&c, i, j);
                    b[[4, j * self.freedom() + 1]] = self.dz(&c, i, j);
                    b[[5, j * self.freedom() + 2]] = self.dx(&c, i, j);
                    b[[4, j * self.freedom() + 2]] = self.dy(&c, i, j);
                    b[[2, j * self.freedom() + 2]] = self.dz(&c, i, j);
                }
                let e = b.dot(u);
                let s = self.d().dot(&e);
                for j in 0..6 {
                    res[[j, i]] += e[j];
                    res[[j + 6, i]] += s[j];
                }
            }
            Ok(res)
        }
    }
    
    impl FiniteElement3D for FE3D4 {
        fn dx(&self, c: &Array2<f64>, _i: usize, j: usize) -> f64 {
            c[[1, j]]
        }
        fn dy(&self, c: &Array2<f64>, _i: usize, j: usize) -> f64 {
            c[[2, j]]
        }
        fn dz(&self, c: &Array2<f64>, _i: usize, j: usize) -> f64 {
            c[[3, j]]
        }
        fn dxi(&self, _i: usize, j: usize) -> f64 {
            [ -1.0, 1.0, 0.0, 0.0 ][j]
        }
        fn deta(&self, _i: usize, j: usize) -> f64 {
            [ -1.0, 0.0, 1.0, 0.0 ][j]
        }
        fn dpsi(&self, _i: usize, j: usize) -> f64 {
            [ -1.0, 0.0, 0.0, 1.0 ][j]
        }
    }
    
    impl FiniteElement3D for FE3D8 {
        fn dx(&self, c: &Array2<f64>, i: usize, j: usize) -> f64 {
            c[[1, j]] + c[[4, j]] * self.x[[i, 1]] + c[[5, j]] * self.x[[i, 2]] + c[[7, j]] * self.x[[i, 1]] * self.x[[i, 2]]
        }
        fn dy(&self, c: &Array2<f64>, i: usize, j: usize) -> f64 {
            c[[2, j]] + c[[4, j]] * self.x[[i, 0]] + c[[6, j]] * self.x[[i, 2]] + c[[7, j]] * self.x[[i, 0]] * self.x[[i, 2]]
        }
        fn dz(&self, c: &Array2<f64>, i: usize, j: usize) -> f64 {
            c[[3, j]] + c[[5, j]] * self.x[[i, 0]] + c[[6, j]] * self.x[[i, 1]] + c[[7, j]] * self.x[[i, 0]] * self.x[[i, 1]]
        }
        fn dxi(&self, i: usize, j: usize) -> f64 {
            let eta = array![ -0.57735026919, -0.57735026919, 0.57735026919, 0.57735026919, -0.57735026919, -0.57735026919, 0.57735026919, 0.57735026919 ];
            let psi = array![ -0.57735026919, 0.57735026919, -0.57735026919, 0.57735026919, -0.57735026919, 0.57735026919, -0.57735026919, 0.57735026919 ];
            [ -0.125 * (1.0 - eta[i]) * (1.0 - psi[i]), 0.125 * (1.0 - eta[i]) * (1.0 - psi[i]),
               0.125 * (1.0 + eta[i]) * (1.0 - psi[i]), -0.125 * (1.0 + eta[i]) * (1.0 - psi[i]),
              -0.125 * (1.0 - eta[i]) * (1.0 + psi[i]), 0.125 * (1.0 - eta[i]) * (1.0 + psi[i]),
               0.125 * (1.0 + eta[i]) * (1.0 + psi[i]), -0.125 * (1.0 + eta[i]) * (1.0 + psi[i]) ][j]
        }
        fn deta(&self, i: usize, j: usize) -> f64 {
            let xi = array![ -0.57735026919, -0.57735026919, -0.57735026919, -0.57735026919, 0.57735026919, 0.57735026919, 0.57735026919, 0.57735026919 ];
            let psi = array![ -0.57735026919, 0.57735026919, -0.57735026919, 0.57735026919, -0.57735026919, 0.57735026919, -0.57735026919, 0.57735026919 ];
            [ -0.125 * (1.0 - xi[i]) * (1.0 - psi[i]), -0.125 * (1.0 + xi[i]) * (1.0 - psi[i]),
               0.125 * (1.0 + xi[i]) * (1.0 - psi[i]),  0.125 * (1.0 - xi[i]) * (1.0 - psi[i]),
              -0.125 * (1.0 - xi[i]) * (1.0 + psi[i]), -0.125 * (1.0 + xi[i]) * (1.0 + psi[i]),
               0.125 * (1.0 + xi[i]) * (1.0 + psi[i]),  0.125 * (1.0 - xi[i]) * (1.0 + psi[i]) ][j]
        }
        fn dpsi(&self, i: usize, j: usize) -> f64 {
            let xi = array![ -0.57735026919, -0.57735026919, -0.57735026919, -0.57735026919, 0.57735026919, 0.57735026919, 0.57735026919, 0.57735026919 ];
            let eta = array![ -0.57735026919, -0.57735026919, 0.57735026919, 0.57735026919, -0.57735026919, -0.57735026919, 0.57735026919, 0.57735026919 ];
            [ -0.125 * (1.0 - xi[i]) * (1.0 - eta[i]), -0.125 * (1.0 + xi[i]) * (1.0 - eta[i]),
              -0.125 * (1.0 + xi[i]) * (1.0 + eta[i]), -0.125 * (1.0 - xi[i]) * (1.0 + eta[i]),
               0.125 * (1.0 - xi[i]) * (1.0 - eta[i]),  0.125 * (1.0 + xi[i]) * (1.0 - eta[i]),
               0.125 * (1.0 + xi[i]) * (1.0 + eta[i]),  0.125 * (1.0 - xi[i]) * (1.0 + eta[i]) ][j]
        }
    }
       
    impl FiniteElement for FE3D4 {
        fn w(&self) -> Array1<f64> {
            array![ -0.13333333333, 0.075, 0.075, 0.075, 0.075 ]
        }
        fn e(&self) -> Array1<f64> {
            array![ self.e[0], self.e[1] ]
        }
        fn x(&self, i: usize, j: usize) -> f64 {
            self.x[[i, j]]
        }
        fn size(&self) -> usize {
            4
        }
        fn freedom(&self) -> usize {
            3
        }
        fn shape_factor(&self, i: usize, j: usize) -> f64 {
            [ 1.0, self.x[[i, 0]], self.x[[i, 1]], self.x[[i, 2]] ][j]
        }
    }
    
    impl FiniteElement for FE3D8 {
        fn w(&self) -> Array1<f64> {
            array![ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ]
        }
        fn e(&self) -> Array1<f64> {
            array![ self.e[0], self.e[1] ]
        }
        fn x(&self, i: usize, j: usize) -> f64 {
            self.x[[i, j]]
        }
        fn size(&self) -> usize {
            8
        }
        fn freedom(&self) -> usize {
            3
        }
        fn shape_factor(&self, i: usize, j: usize) -> f64 {
            [ 1.0, self.x[[i, 0]], self.x[[i, 1]], self.x[[i, 2]], self.x[[i, 0]] * self.x[[i, 1]], 
              self.x[[i, 0]] * self.x[[i, 2]], self.x[[i, 1]] * self.x[[i, 2]], self.x[[i, 0]] * self.x[[i, 1]] * self.x[[i, 2]] ][j]
        }
    }
    
    
    impl FE3D4 {
        pub fn new(x: Array2<f64>, e: [f64; 2]) -> Self {
            Self { x, e, }
        }
    }
    
    impl FE3D8 {
        pub fn new(x: Array2<f64>, e: [f64; 2]) -> Self {
            Self { x, e, }
        }
    }
}