#[path = "core/fem.rs"]
mod fem;
#[path = "core/error.rs"]
mod error;

// 1. Изменить процедуру учета граничных условий (обрабатывать только ненулевые элементы)
// 2. Выборка индексов КЭ (переделать) (x)
// 3. Парсер - предкомпиляция (x)
// 4. Параллельное умножение матрицы на число
// 5. unwrap() (Parser)
// 6. Sparse - попробовать крайт sprs
// 7. Запись в mesh-файл информации о связях
// 8. Вывод на экран информации о погрешности при итерационном решении СЛАУ

//use error::Error;
use fem::Direct;

#[allow(dead_code)]
fn test_1d2() {
    let mesh_name = "D:/Work/python/pyfem/mesh/body1d.trpa";
    // let mesh_name = "/home/serg/work/python/pyfem/mesh/body1d.trpa";
    let mut fem: fem::FEM = match fem::FEM::new(mesh_name) {
        Err(err) => {
            println!("{}", err.say_error());
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_young_modulus(203200.);
    // fem.set_poisons_ratio(0.33);
    fem.set_thickness(1.0);
    fem.set_eps(1.0e-6);
    fem.add_boundary_condition("0", "x == 0", Direct::X);
    fem.add_concentrated_load("1", "x == 10", Direct::X);
    match fem.generate() {
        Err(err) => {
            println!("{}", err.say_error());
            return;
        }
        _ => println!("Done"),
    }
}

#[allow(dead_code)]
fn test_2d4() {
    // let mesh_name = "D:/Work/Qt/QFEM/mesh/console/console_test.trpa";
    // let mesh_name = "D:/Work/Qt/QFEM/mesh/console/console4.trpa";
    // let mesh_name = "D:/Work/python/pyfem/mesh/console.trpa";
    let mesh_name = "/home/serg/work/python/pyfem/mesh/console4.trpa";
    let mut fem: fem::FEM = match fem::FEM::new(mesh_name) {
        Err(err) => {
            println!("{}", err.say_error());
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_young_modulus(203200.0);
    fem.set_poisons_ratio(0.27);
    fem.set_thickness(1.0);
    fem.add_boundary_condition("0", "x == 0", Direct::X | Direct::Y);
    // fem.add_concentrated_load(String::from("-1"), String::from("x == 10"), Direct::Y);
    fem.add_volume_load("1", "", Direct::Y);
    match fem.generate() {
        Err(err) => {
            println!("{}", err.say_error());
            return;
        }
        _ => println!("Done"),
    }
}

#[allow(dead_code)]
fn test_3d4() {
    let mesh_name = "/home/serg/work/Qt/QFEM/mesh/balka.trpa";
    // let mesh_name = "D:/Work/Qt/QFEM/mesh/balka.trpa";
    let mut fem: fem::FEM = match fem::FEM::new(mesh_name) {
        Err(err) => {
            println!("{}", err.say_error());
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_young_modulus(203200.);
    fem.set_poisons_ratio(0.27);
    fem.add_boundary_condition("0", "y == 0", Direct::X | Direct::Y | Direct::Z);
    // fem.add_concentrated_load(String::from("-1"), String::from("y == 4"), Direct::Y);
    fem.add_volume_load("-1", "", Direct::Y);
    match fem.generate() {
        Err(err) => {
            println!("{}", err.say_error());
            return;
        }
        _ => println!("Done"),
    }
}

#[allow(dead_code)]
fn test_3d8() {
    // let mesh_name = "/home/serg/work/mesh/cube_test.trpa";
    // let mesh_name = "D:/Work/python/pyfem/mesh/cube.trpa";
    let mesh_name = "/home/homeniuk/work/python/pyfem/mesh/cube.trpa";
    let mut fem: fem::FEM = match fem::FEM::new(mesh_name) {
        Err(err) => {
            println!("{}", err.say_error());
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_young_modulus(203200.);
    fem.set_poisons_ratio(0.27);
    
    // fem.add_boundary_condition("0", "z == 0", Direct::X | Direct::Y | Direct::Z);
    fem.add_boundary_condition_fun(|_x, _y, _z| { 0. }, |_x, _y, z| { if z == 0. { true } else { false } }, Direct::X | Direct::Y | Direct::Z);
    
    // fem.add_concentrated_load(String::from("1"), String::from("z == 1"), Direct::Z);
    // fem.add_volume_load("1", "", Direct::Z);
    // fem.add_surface_load("1", "z == 1", Direct::Z);
    
    // fem.add_pressure_load("1", "z == 1");

    // fem.add_surface_load_fun(|_x, _y, _z| { -0.5 }, |_x, _y, z| { if z == 1. { true } else { false } }, Direct::Z);
    fem.add_volume_load_fun(|_x, _y, _z| { -0.5 }, |_x, _y, _z| { true }, Direct::Z);
    match fem.generate() {
        Err(err) => {
            println!("{}", err.say_error());
            return;
        }
        _ => println!("Done"),
    }
}

#[allow(dead_code)]
fn test_cyl3d() {
    let mesh_name = "/home/homeniuk/work/python/pyfem/mesh/cyl.trpa";
    let mut fem: fem::FEM = match fem::FEM::new(mesh_name) {
        Err(err) => {
            println!("{}", err.say_error());
            return;
        }
        Ok(fem) => fem,
    };
    fem.set_young_modulus(203200.);
    fem.set_poisons_ratio(0.27);
    
    fem.add_boundary_condition_fun(|_x, _y, _z| { 0. }, |x, _y, _z| { if x == 0.0 { true } else { false } }, Direct::X | Direct::Y | Direct::Z);
    fem.add_boundary_condition_fun(|_x, _y, _z| { 0. }, |x, _y, _z| { if x == 2.0 { true } else { false } }, Direct::X | Direct::Y | Direct::Z);
    
    fem.add_pressure_load_fun(|_x, _y, _z| { 1.0E+4 }, |_x, y, z| { if (y * y + z * z - 0.5 * 0.5).abs() < 1.0E-2 {true} else {false} });
    match fem.generate() {
        Err(err) => {
            println!("{}", err.say_error());
            return;
        }
        _ => println!("Done"),
    }
}


fn main() {
    // test_1d2();
    // test_2d4();
    // test_3d4();
    test_3d8();
    // test_cyl3d();
}
