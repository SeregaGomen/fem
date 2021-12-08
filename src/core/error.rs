// Типы ошибок
#[derive(Debug)]
#[allow(dead_code)]
pub enum Error {
    OpenFile,
    ReadFile,
    WriteFile,
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
    pub fn say_error(&self) -> &str {
        match self {
            Error::OpenFile => "Error: unable open file",
            Error::ReadFile => "Error: unable read file",
            Error::WriteFile => "Error: unable write file",
            Error::InvalidFEType => "Error: invalid FE type",
            Error::InvalidNumber => "Error: invalid number",
            Error::SingularMatrix => "Error: singular matrix",
            Error::InverseMatrix => "Error calculation of inverse matrix",
            Error::DeterminantMatrix => "Error calculation of matrix determinant",
            Error::InvalidIndex => "Error: invalid index",
            Error::UndefError => "Error: undefined variable",
            Error::BracketError => "Error: unbalanced brackets",
            Error::SyntaxError => "Syntax error",
            Error::InternalError => "Internal error",
        }
    }    
}