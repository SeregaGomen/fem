pub mod fem;
pub mod json;


// 1. Изменить процедуру учета граничных условий (обрабатывать только ненулевые элементы)
// 2. Выборка индексов КЭ (переделать) (x)
// 3. Парсер - предкомпиляция (x)
// 4. Параллельное умножение матрицы на число (x)
// 5. unwrap() (Parser)
// 6. Sparse - попробовать крайт sprs
// 7. Запись в mesh-файл информации о связях (x)
// 8. Вывод на экран информации о погрешности при итерационном решении СЛАУ
// 9. Запись результатов (x)