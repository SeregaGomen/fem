//use std;
use std::fmt;
// use std::error::Error;
use std::convert::From;
use std::io::Error as IoError;
use std::num::ParseIntError as ParseIntError;
use std::num::ParseFloatError as ParseFloatError;
use russell_sparse::StrError;
use json::Error as JsonError;

#[derive(Debug)]
pub enum FemError {
    Io(IoError),
    ParseInt(ParseIntError),
    ParseFloat(ParseFloatError),
    JsonError(JsonError),
    StrError(StrError),
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
    YoungModulusError,
    PoissonRatioError,
    ResultError,
    DirectError,
    ValueError,
    PredicateError,
    IncorrectDirectError,
    IncorrectParamError,
    StressStrainCurveError,
    ThicknessError,
    Other,
}

impl fmt::Display for FemError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            FemError::Io(ref cause) => write!(f, "I/O error: {}", cause),
            FemError::ParseInt(ref cause) => write!(f, "Parse int error: {}", cause),
            FemError::ParseFloat(ref cause) => write!(f, "Parse float error: {}", cause),
            FemError::JsonError(ref cause) => write!(f, "Json error: {}", cause),
            FemError::StrError(ref cause) => write!(f, "Russell error: {}", *cause),
            FemError::OpenFile => write!(f, "Unable open file"),
            FemError::ReadFile => write!(f, "Unable read file"),
            FemError::WriteFile => write!(f, "Unable write file"),
            FemError::InvalidFEType => write!(f, "Invalid FE type"),
            FemError::InvalidNumber => write!(f, "Invalid number"),
            FemError::SingularMatrix => write!(f, "Singular matrix"),
            FemError::InverseMatrix => write!(f, "Error calculation of inverse matrix"),
            FemError::DeterminantMatrix => write!(f, "Error calculation of matrix determinant"),
            FemError::InvalidIndex => write!(f, "Invalid index"),
            FemError::UndefError => write!(f, "Undefined variable"),
            FemError::BracketError => write!(f, "Unbalanced brackets"),
            FemError::SyntaxError => write!(f, "Syntax error"),
            FemError::MeshError => write!(f, "Mesh-file not specified"),
            FemError::ResultError => write!(f, "Result-file not specified"),
            FemError::ParamError => write!(f, "Wrong parameter"), 
            FemError::InternalError => write!(f, "Internal error"),
            FemError::YoungModulusError => write!(f, "Modulus of elasticity not set"),
            FemError::PoissonRatioError => write!(f, "Poisson's ratio not set"),
            FemError::DirectError => write!(f, "Direct not set"),
            FemError::IncorrectDirectError => write!(f, "Incorrect direct"),
            FemError::IncorrectParamError => write!(f, "Incorrect parameter"),
            FemError::ValueError => write!(f, "Value not set"),
            FemError::PredicateError => write!(f, "Predicate not set"),
            FemError::StressStrainCurveError => write!(f, "Stress-strain curve not set"),
            FemError::ThicknessError => write!(f, "FE thickness not set"),
            FemError::Other => write!(f, "Unknown error"),
        }
    }
}

// impl Error for FemError {
//     fn cause(&self) -> Option<&dyn Error> {
//         match *self {
//             FemError::Io(ref cause) => Some(cause),
//             FemError::ParseInt(ref cause) => Some(cause),
//             FemError::ParseFloat(ref cause) => Some(cause),
//             FemError::JsonError(ref cause) => Some(cause),
//             FemError::StrError(ref _cause) => Some(&FemError::RussellError), //Some(cause),
//             FemError::OpenFile => Some(&FemError::OpenFile),
//             FemError::ReadFile => Some(&FemError::ReadFile),
//             FemError::WriteFile => Some(&FemError::WriteFile),
//             FemError::InvalidFEType => Some(&FemError::InvalidFEType),
//             FemError::InvalidNumber => Some(&FemError::InvalidNumber),
//             FemError::SingularMatrix => Some(&FemError::SingularMatrix),
//             FemError::InverseMatrix => Some(&FemError::InverseMatrix),
//             FemError::DeterminantMatrix => Some(&FemError::DeterminantMatrix),
//             FemError::InvalidIndex => Some(&FemError::InvalidIndex),
//             FemError::UndefError => Some(&FemError::UndefError),
//             FemError::BracketError => Some(&FemError::BracketError),
//             FemError::SyntaxError => Some(&FemError::SyntaxError),
//             FemError::MeshError => Some(&FemError::MeshError),
//             FemError::ResultError => Some(&FemError::ResultError),
//             FemError::ParamError => Some(&FemError::ParamError),
//             FemError::InternalError => Some(&FemError::InternalError),
//             FemError::YoungModulusError => Some(&FemError::YoungModulusError),
//             FemError::PoissonRatioError => Some(&FemError::PoissonRatioError),
//             FemError::DirectError => Some(&FemError::DirectError),
//             FemError::IncorrectDirectError => Some(&FemError::IncorrectDirectError),
//             FemError::IncorrectParamError => Some(&FemError::IncorrectParamError),
//             FemError::ValueError => Some(&FemError::ValueError),
//             FemError::PredicateError => Some(&FemError::PredicateError),
//             FemError::StressStrainCurveError => Some(&FemError::StressStrainCurveError),
//             FemError::ThicknessError => Some(&FemError::ThicknessError),
//             FemError::RussellError => Some(&FemError::RussellError),
//             FemError::Other => None,
//         }
//     }
// }

// Support converting system errors into our custom error.
// This trait is used in `try!`.
impl From<IoError> for FemError {
    fn from(cause: IoError) -> FemError {
        FemError::Io(cause)
    }
}

impl From<ParseIntError> for FemError {
    fn from(cause: ParseIntError) -> FemError {
        FemError::ParseInt(cause)
    }
}

impl From<ParseFloatError> for FemError {
    fn from(cause: ParseFloatError) -> FemError {
        FemError::ParseFloat(cause)
    }
}

impl From<JsonError> for FemError {
    fn from(cause: JsonError) -> FemError {
        FemError::JsonError(cause)
    }
}

impl From<StrError> for FemError {
    fn from(cause: StrError) -> FemError {
        FemError::StrError(cause)
    }
}
