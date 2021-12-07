// Типы ошибок
#[derive(Debug)]
#[allow(dead_code)]
pub enum Error {
    OpenFile,
    ReadFile,
    InvalidFEType, 
    InvalidNumber,
    SingularMatrix,  
    InverseMatrix,  
    DeterminantMatrix,
    InvalidIndex,
    UndefError,
    BracketError,
    SyntaxError,
    InternalError,
}
#[allow(dead_code)]
impl Error {
    pub fn say_error(&self) -> String {
        match self {
            Error::OpenFile => String::from("Error: unable open mesh-file"),
            Error::ReadFile => String::from("Error: unable read mesh-file"),
            Error::InvalidFEType => String::from("Error: invalid FE type"),
            Error::InvalidNumber => String::from("Error: invalid number"),
            Error::SingularMatrix => String::from("Error: singular matrix"),
            Error::InverseMatrix => String::from("Error calculation of inverse matrix"),
            Error::DeterminantMatrix => String::from("Error calculation of matrix determinant"),
            Error::InvalidIndex => String::from("Error: invalid index"),
            Error::UndefError => String::from("Error: undefined variable"),
            Error::BracketError => String::from("Error: unbalanced brackets"),
            Error::SyntaxError => String::from("Syntax error"),
            Error::InternalError => String::from("Internal error"),
        }
    }    
}