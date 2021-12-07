// Утилиты
use ndarray::prelude::*;
use crate::error::Error;

#[allow(dead_code)]
// Решение СЛАУ методом Гаусса
pub fn solve(a: &mut Array2<f64>, b: &mut Array1<f64>, eps: f64) -> Result<Array1<f64>, Error> {
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
        return Err(Error::SingularMatrix);   
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
pub fn inv(m: &Array2<f64>) -> Result<Array2<f64>, Error> {
    let ret;
    if m.shape()[0] == 1 && m.shape()[1] == 1 {
        ret = array![[1.0 / m[[0, 0]]]];
    }
    else if m.shape()[0] == 2 && m.shape()[1] == 2 {
        ret = array![[m[[1, 1]], -m[[0, 1]]], [-m[[1, 0]], m[[0, 0]]]] / det(m)?;
    }
    else if m.shape()[0] == 3 && m.shape()[1] == 3 {
        ret = array![[m[[1, 1]] * m[[2, 2]] - m[[1, 2]] * m[[2, 1]], m[[0, 2]] * m[[2, 1]] - m[[0, 1]] * m[[2, 2]], m[[0, 1]] * m[[1, 2]] - m[[0, 2]] * m[[1, 1]]], 
                    [m[[1, 2]] * m[[2, 0]] - m[[1, 0]] * m[[2, 2]], m[[0, 0]] * m[[2, 2]] - m[[0, 2]] * m[[2, 0]], m[[0, 2]] * m[[1, 0]] - m[[0, 0]] * m[[1, 2]]], 
                    [m[[1, 0]] * m[[2, 1]] - m[[1, 1]] * m[[2, 0]], m[[0, 1]] * m[[2, 0]] - m[[0, 0]] * m[[2, 1]], m[[0, 0]] * m[[1, 1]] - m[[0, 1]] * m[[1, 0]]]] / det(m)?;
    }
    else {
        return Err(Error::InverseMatrix);
    }
    Ok(ret)
}

// Определитель матрицы
pub fn det(m: &Array2<f64>) -> Result<f64, Error> {
    let ret;
    if m.shape()[0] == 1 && m.shape()[1] == 1 {
        ret = m[[0, 0]];
    }
    else if m.shape()[0] == 2 && m.shape()[1] == 2 {
        ret = m[[0, 0]] * m[[1, 1]] - m[[0, 1]] * m[[1, 0]];
    }
    else if m.shape()[0] == 3 && m.shape()[1] == 3 {
        ret = m[[0, 0]] * m[[1, 1]] * m[[2, 2]] + m[[0, 1]] * m[[1, 2]] * m[[2, 0]] + m[[0, 2]] * m[[1, 0]] * m[[2, 1]] -
                m[[0, 2]] * m[[1, 1]] * m[[2, 0]] - m[[0, 0]] * m[[1, 2]] * m[[2, 1]] - m[[0, 1]] * m[[1, 0]] * m[[2, 2]];
    }
    else {
        return Err(Error::DeterminantMatrix);
    }
    Ok(ret)
}

// Скалярное произведение векторов
pub fn scalar_product(lhs: &Array1<f64>, rhs: &Array1<f64>) -> f64 {
    let mut ret = 0.;
    for i in 0..lhs.len() {
        ret += lhs[i] * rhs[i];
    }
    ret
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