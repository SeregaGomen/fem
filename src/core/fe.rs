// Поддержка конечного элемента
use std::fmt;
use ndarray::prelude::*;
use crate::fem::util;
use crate::error::Error;

// Типы конечных элементов (КЭ)
#[derive(Copy, Clone, PartialEq)]
pub enum FEType {
    FE1D2,
    FE2D3,
    FE2D4,
    FE3D4,
    FE3D8,
}

impl fmt::Display for FEType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let str = match self {
            FEType::FE1D2 => "fe1d2",     
            FEType::FE2D3 => "fe2d3",     
            FEType::FE2D4 => "fe2d4",     
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
    fn create(&self) -> Result<Array2<f64>, Error> {
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
    use crate::error::Error; 
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
        fn generate(&mut self) -> Result<Array2<f64>, Error> {
            // self.create()?;
            let size = self.size() * self.freedom();
            let mut local: Array2<f64> = Array2::zeros((size , size));
            let mut b: Array2<f64> = Array2::zeros((1, self.size()));
            for i in 0..self.w().len() {
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
        fn calc(&self, u: &Array1<f64>) -> Result<Array2<f64>, Error> {
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
        pub fn new(e: [f64; 2], thk: f64, x: Array2<f64>) -> Self {
            Self { e, thk, x, }
        }
    }
}

pub mod fe2d {
    use crate::error::Error; 
    use super::FiniteElement;
    use ndarray::prelude::*;
    use crate::fem::util;

    pub struct FE2D3 {
        e: [f64; 2],        
        thk: f64,           
        x: Array2<f64>,     
    }
    
    pub struct FE2D4 {
        e: [f64; 2],        
        thk: f64,           
        x: Array2<f64>,     
    }
    
    pub trait FiniteElement2D: FiniteElement {
        fn thk(&self) -> f64;
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
        fn generate(&mut self) -> Result<Array2<f64>, Error> {
            // self.create()?;
            let size = self.size() * self.freedom();
            let mut local: Array2<f64> = Array2::zeros((size , size));
            let mut b: Array2<f64> = Array2::zeros((3, self.size() * self.freedom()));
            for i in 0..self.w().len() {
                // Якобиан и обратная матрица Якоби
                let jacobi = self.jacobi(i);
                let inv_jacobi = util::inv(&jacobi)?;
                let jacobian = util::det(&jacobi)?;
                for j in 0..self.size() {
                    b[[0, j * self.freedom() + 0]] = inv_jacobi[[0, 0]] * self.dxi(i, j) + inv_jacobi[[0, 1]] * self.deta(i, j);
                    b[[1, j * self.freedom() + 1]] = inv_jacobi[[1, 0]] * self.dxi(i, j) + inv_jacobi[[1, 1]] * self.deta(i, j);
                    b[[2, j * self.freedom() + 1]] = b[[0, j * self.freedom() + 0]];
                    b[[2, j * self.freedom() + 0]] = b[[1, j * self.freedom() + 1]];
                }
                local = local + b.t().dot(&self.d()).dot(&b) * (self.w()[i] * self.thk() * jacobian.abs());
            }
            Ok(local)
        }
        fn calc(&self, u: &Array1<f64>) -> Result<Array2<f64>, Error> {
            let c = self.create()?;
            let mut res: Array2<f64> = Array2::zeros((6, self.size() * self.freedom()));
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
            2
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
            2
        }
        fn shape_factor(&self, i: usize, j: usize) -> f64 {
            [ 1.0, self.x[[i, 0]], self.x[[i, 1]] ][j]
        }
    }
    
    impl FiniteElement2D for FE2D3 {
        fn thk(&self) -> f64 {
            self.thk
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
    }
    
    impl FiniteElement2D for FE2D4 {
        fn thk(&self) -> f64 {
            self.thk
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
    }

    impl FE2D3 {
        pub fn new(e: [f64; 2], thk: f64, x: Array2<f64>) -> Self {
            Self { e, thk, x, }
        }
    }
    
    impl FE2D4 {
        pub fn new(e: [f64; 2], thk: f64, x: Array2<f64>) -> Self {
            Self { e, thk, x, }
        }
    }
}

pub mod fe3d {
    use crate::error::Error; 
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
        fn generate(&mut self) -> Result<Array2<f64>, Error> {
            // self.create()?;
            let size = self.size() * self.freedom();
            let mut local: Array2<f64> = Array2::zeros((size , size));
            let mut b: Array2<f64> = Array2::zeros((6, self.size() * self.freedom()));
            for i in 0..self.w().len() {
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
        fn calc(&self, u: &Array1<f64>) -> Result<Array2<f64>, Error> {
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
        pub fn new(e: [f64; 2], x: Array2<f64>) -> Self {
            Self { e, x, }
        }
    }
    
    impl FE3D8 {
        pub fn new(e: [f64; 2], x: Array2<f64>) -> Self {
            Self { e, x, }
        }
    }
}





