//use std;
use std::fmt;

// Типы ошибок
#[derive(Debug)]
#[allow(dead_code)]
pub enum ErrorCode {
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
    MeshError,
    ParamError,
    JsonError,
    YoungModulusError,
    PoissonRatioError,
    ResultError,
    DirectError,
    ValueError,
    PredicateError,
    IncorrectDirectError,
}

#[derive(Debug)]
pub struct Error { 
    code: ErrorCode, 
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> { 
        write!(f, "{}", self.message())
    }
}

impl std::error::Error for Error {
    fn description(&self) -> &str {
        self.message()
    }
}

#[allow(dead_code)]
impl Error {
    pub fn new(code: ErrorCode) -> Self {
        Self{code}
    }
    fn message(&self) -> &str {
        match self.code {
            ErrorCode::OpenFile => "Unable open file",
            ErrorCode::ReadFile => "Unable read file",
            ErrorCode::WriteFile => "Unable write file",
            ErrorCode::InvalidFEType => "Invalid FE type",
            ErrorCode::InvalidNumber => "Invalid number",
            ErrorCode::SingularMatrix => "Singular matrix",
            ErrorCode::InverseMatrix => "Error calculation of inverse matrix",
            ErrorCode::DeterminantMatrix => "Error calculation of matrix determinant",
            ErrorCode::InvalidIndex => "Invalid index",
            ErrorCode::UndefError => "Undefined variable",
            ErrorCode::BracketError => "Unbalanced brackets",
            ErrorCode::SyntaxError => "Syntax error",
            ErrorCode::MeshError => "Mesh-file not specified",
            ErrorCode::ResultError => "Result-file not specified",
            ErrorCode::JsonError => "Json-file read error", 
            ErrorCode::ParamError => "Wrong parameter", 
            ErrorCode::InternalError => "Internal error",
            ErrorCode::YoungModulusError => "Modulus of elasticity not set",
            ErrorCode::PoissonRatioError => "Poisson's ratio not set",
            ErrorCode::DirectError => "Direct not set",
            ErrorCode::IncorrectDirectError => "Incorrect direct",
            ErrorCode::ValueError => "Value not set",
            ErrorCode::PredicateError => "Predicate not set",
        }
    }
}

#[allow(dead_code)]
pub fn error(code: ErrorCode) -> Error {
    Error::new(code)
}