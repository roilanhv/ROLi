% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funcoes_roli.R
\name{correct.clearness.index}
\alias{correct.clearness.index}
\title{Create from hourly clearness index a mean daily clearness index and 
smooth clearness index}
\usage{
correct.clearness.index(data_, lon, lat, timezone, window_size = 6)
}
\arguments{
\item{data_}{Dataframe with at least date and Rg column.}

\item{lon}{longitude from local analisys}

\item{lat}{latitude from local analisys}

\item{timezone}{local time diference with GMT (-1 for fluxes measurement)}

\item{window_size}{lenght of window to smoothing}
}
\value{
The same input dataframe with K_hourly, K_mean (daily mean), 
K_smooth (smooth K_mean index in specific temporal window) with K index columns
}
\description{
Create from hourly clearness index a mean daily clearness index and 
smooth clearness index
}
\author{
Roilan Hernandez
}
