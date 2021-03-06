% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions2eval.R
\name{get.AllSchems}
\alias{get.AllSchems}
\title{Function to get all scheme calculation, adjust = TRUE make a adjusting NLS.}
\usage{
get.AllSchems(data, Ovrcst_sch = c("CQB", "CKC", "CCB", "CKZ", "CWU", "CJG",
  "CLM", "CFG"), Emiss_sch = c("EAN", "EBR", "EBT", "EDO", "EGR", "EIJ",
  "EID", "EKZ", "ENM", "EPR", "EST", "ESW", "EAI"), adjust = FALSE,
  method = "non-linear", nsample = 1000, max_iter = 10, stats = "rmse",
  log_file = NULL, forced = TRUE)
}
\arguments{
\item{data}{Data frame with atmospheric variables}

\item{Ovrcst_sch}{Schemes for cloud cover index}

\item{Emiss_sch}{Schemes of atmosphere emissivity}

\item{adjust}{FALSE, TRUE for NLS adjusting}

\item{method}{"non-linear" (default) for Non linear Least Square adjust,
"montecarlo" for MonteCarlo optimization. Later be usseful when a NLS can't 
adjust observed data allowing optimization.}

\item{nsample}{population number evaluated in each iteration 
(only when method = "montecarlo").}

\item{max_iter}{maximun number of iterations (only when method = "montecarlo").}

\item{stats}{statistical function to be minimized (only when method = "montecarlo"),
NOTE: the best result should be 0.0 (ex., if stats = r (correlation), then transform to 
rMod = 1.0 - r, so the best result is when r== 1.0, so rMod == 0.0)}

\item{log_file}{path to a log file}

\item{forced}{Forced adjust, default TRUE}
}
\value{
Data frame with Li observed and all combintions of schemes for calculations of Li
}
\description{
Function to get all scheme calculation, adjust = TRUE make a adjusting NLS.
}
\author{
Roilan Hernandez
}
