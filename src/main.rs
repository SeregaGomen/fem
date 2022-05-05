#[path = "fem/fem.rs"]
mod fem;
#[path = "fem/error.rs"]
mod error;
#[path = "fem/json.rs"]
mod json;


// 1. Изменить процедуру учета граничных условий (обрабатывать только ненулевые элементы)
// 2. Выборка индексов КЭ (переделать) (x)
// 3. Парсер - предкомпиляция (x)
// 4. Параллельное умножение матрицы на число (x)
// 5. unwrap() (Parser)
// 6. Sparse - попробовать крайт sprs
// 7. Запись в mesh-файл информации о связях (x)
// 8. Вывод на экран информации о погрешности при итерационном решении СЛАУ
// 9. Запись результатов (x)

use fem::Direct;

#[allow(dead_code)]
fn test_1d2(nthreads: usize) {
    let file_name = ("data/body1d.mesh", "data/body1d.res");
    let mut fem: fem::FEM = match fem::FEM::new(file_name.0) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_num_threads(nthreads);
    fem.set_young_modulus(203200.);
    fem.set_thickness(1.0);
    fem.set_eps(1.0e-6);
    fem.add_boundary_condition("0", "x == 0", Direct::X);
    fem.add_concentrated_load("1", "x == 10", Direct::X);
    match fem.generate(file_name.1) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        _ => println!("Done"),
    }
}

#[allow(dead_code)]
fn test_2d4(nthreads: usize) {
    let file_name = ("data/console.mesh", "data/console.res");
    let mut fem: fem::FEM = match fem::FEM::new(file_name.0) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_num_threads(nthreads);
    fem.set_young_modulus(203200.0);
    fem.set_poisons_ratio(0.27);

    
    fem.set_thickness(1.0);
    fem.add_boundary_condition("0", "x == 0", Direct::X | Direct::Y);
    fem.add_concentrated_load("-1", "x == 10", Direct::Y);
    match fem.generate(file_name.1) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        _ => println!("Done"),
    }
}

#[allow(dead_code)]
fn test_3d8(nthreads: usize) {
    let file_name = ("data/cube.mesh", "data/cube.res");
    let mut fem: fem::FEM = match fem::FEM::new(file_name.0) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_num_threads(nthreads);
    fem.set_young_modulus(203200.);
    fem.set_poisons_ratio(0.27);
    
    // fem.add_boundary_condition("0", "z == 0", Direct::X | Direct::Y | Direct::Z);
    fem.add_boundary_condition_fun(|_x, _y, _z| { 0. }, |_x, _y, z| { if z == 0. { true } else { false } }, Direct::X | Direct::Y | Direct::Z);
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
fn test_shell_3(nthreads: usize) {
    let file_name = ("data/shell-tube-3.mesh", "data/shell-tube-3.res");
    let mut fem: fem::FEM = match fem::FEM::new(file_name.0) {
        Err(err) => {
            println!("{}", err);
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_num_threads(nthreads);
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
    test_1d2(8);
    // test_2d4(8);
    // test_3d8(8);
    // test_shell_3(8);
}


// fn main() {
//     use std::env;

//     let args: Vec<String> = env::args().collect();
    
//     if args.len() < 2 {
//         println!("Too few parameters!");
//         return;
//     }
//     match json::read_json(args[1].as_str()) {
//          Ok(_) => println!("Done"),
//          Err(e) => println!("Error: {}", e),
//     };
// }

