// Утилиты
use std::io::prelude::*;
use std::str;
use ndarray::prelude::*;
use super::error::FemError;


#[allow(dead_code)]
// Решение СЛАУ методом Гаусса
pub fn solve(a: &mut Array2<f64>, b: &mut Array1<f64>, eps: f64) -> Result<Array1<f64>, FemError> {
    let n = a.shape()[0];
    for i in 0..n - 1 {
        if a[[i, i]].abs() < eps {
            for j in i + 1..n {
                if a[[j, i]].abs() < eps {
                    continue;
                }
                for k in 0..n {
                    a.swap([j, k], [i, k]);
                }
                b.swap([j], [i]);
            }
        }
        let k1 = a[[i, i]];
        for j in i + 1..n {
            if a[[j, i]].abs() < eps {
                continue;
            }
            let k2 = a[[j, i]];
            for k in 0..n {
                a[[j, k]] -= k2 * a[[i, k]] / k1;
            }
            b[j] -= k2 * b[i] / k1;
        }
    }
    if a[[n - 1, n - 1]].abs() < eps {
        return Err(FemError::SingularMatrix);   
    }
    let mut x: Array1<f64> = Array1::zeros(n);
    x[n - 1] = b[n - 1] / a[[n - 1, n - 1]];
    for i in (0..n - 1).rev() {
        x[i] = b[i];
        for k in i + 1..n {
            x[i] -= x[k] * a[[i, k]];
        }
        x[i] /= a[[i, i]];
    }
    Ok(x)
}

// Обратная матрица
pub fn inv(m: &Array2<f64>) -> Result<Array2<f64>, FemError> {
    let ret;
    if m.shape()[0] == 1 && m.shape()[1] == 1 {
        ret = array![[1.0 / m[[0, 0]]]];
    } else if m.shape()[0] == 2 && m.shape()[1] == 2 {
        ret = array![[m[[1, 1]], -m[[0, 1]]], [-m[[1, 0]], m[[0, 0]]]] / det(m)?;
    } else if m.shape()[0] == 3 && m.shape()[1] == 3 {
        ret = array![[m[[1, 1]] * m[[2, 2]] - m[[1, 2]] * m[[2, 1]], m[[0, 2]] * m[[2, 1]] - m[[0, 1]] * m[[2, 2]], m[[0, 1]] * m[[1, 2]] - m[[0, 2]] * m[[1, 1]]], 
                    [m[[1, 2]] * m[[2, 0]] - m[[1, 0]] * m[[2, 2]], m[[0, 0]] * m[[2, 2]] - m[[0, 2]] * m[[2, 0]], m[[0, 2]] * m[[1, 0]] - m[[0, 0]] * m[[1, 2]]], 
                    [m[[1, 0]] * m[[2, 1]] - m[[1, 1]] * m[[2, 0]], m[[0, 1]] * m[[2, 0]] - m[[0, 0]] * m[[2, 1]], m[[0, 0]] * m[[1, 1]] - m[[0, 1]] * m[[1, 0]]]] / det(m)?;
    } else {
        return Err(FemError::InverseMatrix);
    }
    Ok(ret)
}

// Определитель матрицы
pub fn det(m: &Array2<f64>) -> Result<f64, FemError> {
    let ret;
    if m.shape()[0] == 1 && m.shape()[1] == 1 {
        ret = m[[0, 0]];
    } else if m.shape()[0] == 2 && m.shape()[1] == 2 {
        ret = m[[0, 0]] * m[[1, 1]] - m[[0, 1]] * m[[1, 0]];
    } else if m.shape()[0] == 3 && m.shape()[1] == 3 {
        ret = m[[0, 0]] * m[[1, 1]] * m[[2, 2]] + m[[0, 1]] * m[[1, 2]] * m[[2, 0]] + m[[0, 2]] * m[[1, 0]] * m[[2, 1]] -
                m[[0, 2]] * m[[1, 1]] * m[[2, 0]] - m[[0, 0]] * m[[1, 2]] * m[[2, 1]] - m[[0, 1]] * m[[1, 0]] * m[[2, 2]];
    } else {
        return Err(FemError::DeterminantMatrix);
    }
    Ok(ret)
}

// Форматированный вывод вещественных чисел
pub fn fmt_f64(num: f64, width: usize, precision: usize, exp_pad: usize) -> String {
    let mut num = format!("{:.precision$e}", num, precision = precision);
    // Safe to `unwrap` as `num` is guaranteed to contain `'e'`
    let exp = num.split_off(num.find('e').unwrap());

    let (sign, exp) = if exp.starts_with("e-") {
        ('-', &exp[2..])
    } else {
        ('+', &exp[1..])
    };
    num.push_str(&format!("e{}{:0>pad$}", sign, exp, pad = exp_pad));
    format!("{:+>width$}", num, width = width)
}

pub fn get_max(row: ArrayView1<f64>) ->f64 {
    let mut max = row[0];
    for i in 1..row.len() {
        if max < row[i] {
            max = row[i];
        }
    }
    max
}

pub fn get_min(row: ArrayView1<f64>) ->f64 {
    let mut min = row[0];
    for i in 1..row.len() {
        if min > row[i] {
            min = row[i];
        }
    }
    min
}

fn norm3(v: &Array1<f64>) -> Array1<f64> {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt(); 
    array![v[0] / len, v[1] / len, v[2] / len]
}

fn create_vector(p1: ArrayView1<f64>, p2: ArrayView1<f64>) -> Array1<f64> {
    let mut res = Array1::<f64>::zeros(3);
    for i in 0..3 {
        res[i] = p2[i] - p1[i];
    }
    norm3(&res)
}

fn cross_product(a: &Array1<f64>, b: &Array1<f64>) -> Array1<f64> {
    let mut res = Array1::<f64>::zeros(3);
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
    norm3(&res)
}

// Построение матрицы преобразования для оболочечных КЭ
pub fn create_transform_matrix(x: &Array2<f64>) -> Array2<f64> {
    let tmp = create_vector(x.row(2), x.row(0));
    let vx = create_vector(x.row(1), x.row(0));
    let vz = cross_product(&vx, &tmp);
    let vy = cross_product(&vz, &vx);
    array![[vx[0], vx[1], vx[2]], [vy[0], vy[1], vy[2]], [vz[0], vz[1], vz[2]]]
}

// Вспомогательная матрица преобразования для оболочечных КЭ
pub fn create_ext_transform_matrix(v: &Array2<f64>, size: usize, freedom: usize) -> Array2<f64> {
    let mut m = Array2::<f64>::zeros((size * freedom, size * freedom));
    for i in 0..3 {
        for j in 0..3 {
            for k in (0..size * freedom).step_by(3) {
                m[[i + k, j + k]] = v[[i, j]];
                m[[i + k, j + k]] = v[[i, j]];
                m[[i + k, j + k]] = v[[i, j]];
            }
        }
    }    
    m
}

// Чтение первой непустой строки из файла
pub fn get_line<R: BufRead>(reader: &mut R) -> Result<String, FemError> {
    let mut res;
    let mut bytes: Vec<u8> = Vec::new();
    loop {
        let len = match reader.read_until(b'\n', &mut bytes) {
            Err(_) => return Err(FemError::ReadFile),
            Ok(len) => len,
        };
        if len == 0 {
            return Ok(String::new());
        }
        res = match str::from_utf8(&bytes) {
            Ok(res) => res,
            Err(_) => return Err(FemError::InvalidNumber),
        };
        if res.trim().len() > 0 {
            break;
        }
    }
    Ok(res.to_string())
}
