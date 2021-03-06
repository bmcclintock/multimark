% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multimark_functions.R
\docType{class}
\name{multimarksetup-class}
\alias{multimarksetup-class}
\title{Class \code{"multimarksetup"}}
\description{
A class of 'mulitmark' model inputs
}
\section{Slots}{

\describe{
\item{\code{Enc.Mat}}{Object of class \code{"matrix"}. The observed encounter histories (with rows corresponding to individuals and columns corresponding to sampling occasions).}

\item{\code{data.type}}{Object of class \code{"character"}. The encounter history data type ("never", "sometimes", or "always").}

\item{\code{vAll.hists}}{Object of class \code{"integer"}. An ordered vector containing all possible encounter histories in sequence.}

\item{\code{Aprime}}{Object of class \code{"sparseMatrix"}. Transpose of the A matrix mapping latent encounter histories to observed histories.}

\item{\code{indBasis}}{Object of class \code{"numeric"}.An ordered vector of the indices of the three encounter histories updated by each basis vector.}

\item{\code{ncolbasis}}{Object of class \code{"integer"}. The number of needed basis vectors.}

\item{\code{knownx}}{Object of class \code{"integer"}. Frequencies of known encounter histories.}

\item{\code{C}}{Object of class \code{"integer"}. Sampling occasion of first capture for each encounter history.}

\item{\code{L}}{Object of class \code{"integer"}. Sampling occasion of last capture for each encounter history.}

\item{\code{naivex}}{Object of class \code{"integer"}. ``Naive'' latent history frequencies assuming a one-to-one mapping with \code{Enc.Mat}.}

\item{\code{covs}}{Object of class \code{"data.frame"}. Temporal covariates for detection probability (the number of rows in the data frame must equal the number of sampling occasions).}
}}

\section{Objects from the Class}{

Objects can be created by calls of the form \code{processdata(Enc.Mat, ...)} or \code{new("multimarksetup", ...)}.
}

\section{Methods}{

No methods defined with class "multimarksetup".
}

\examples{
showClass("multimarksetup")
}
\seealso{
\code{\link{processdata}}
}
\author{
Brett T. McClintock
}
\keyword{classes}
