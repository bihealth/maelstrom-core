use std::{
    num::ParseIntError,
    process::{ExitCode, Termination},
};

#[derive(thiserror::Error, Debug, Clone)]
pub enum AppError {
    // #[error("Internal error.")]
    // Internal,
}

impl Termination for AppError {
    fn report(self) -> ExitCode {
        match self {
            // Internal => ExitCode::from(1),
            //   Other => ExitCode::from(255),
        }
    }
}

#[derive(thiserror::Error, Debug, Clone)]
pub enum ArgError {
    #[error("Invalid format in interval")]
    IntervalInvalidFormat,
    #[error("Invalid integer coordinates in interval")]
    IntervalInvalidInts(#[from] ParseIntError),
}
