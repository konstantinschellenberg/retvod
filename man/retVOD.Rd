% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retVOD.R
\name{retVOD}
\alias{retVOD}
\title{Retrieve VOD for all values}
\usage{
retVOD(
  tbH,
  tbV,
  smc,
  vod,
  Tair,
  Tsoil,
  omega,
  roughness,
  cf,
  inc_angle,
  silent = F
)
}
\arguments{
\item{tbH}{Horizontal brightness temperatures}

\item{tbV}{Vertical brightness temperatures}

\item{smc}{Soil moisture}

\item{vod}{VOD estimate}

\item{Tair}{Air temperature}

\item{Tsoil}{Soil temperature}

\item{omega}{Canopy opacity, scattering albedo}

\item{roughness}{Soil roughness estimate}

\item{cf}{Clay fraction}

\item{inc_angle}{Incidence angle}

\item{silent}{Silence progress bars in the console. Default is FALSE.}
}
\value{
Retrieved VOD and auxillary information
}
\description{
Retrieve VOD for all values
}
