\name{Data}
\alias{Es}
\alias{H2O}
\alias{Parea}
\alias{Parea.merra}
\alias{Parea.merra.glob}
\alias{Parea.month}
\alias{Parea.merra.month}
\alias{Parea.eraint.month}
\alias{Parea.cmip5}
\alias{Qs}
\alias{Rho}
\alias{prtot.merra}
\alias{prtot.eraint}
\alias{raingaugedata.USA}
\alias{raingaugedata.World}
\title{Global data extract from ERAINT and TRMM}
\description{Indices describing the vapour saturation pressure based on
  ERAINT 2-m temperature (\code{Es}), the total air column water content
  from ERAINT (\code{H2O}, the precipitation area from TRMM
  (\code{Parea}), surface moisture flux from ERAINT (\code{Qs}), and the
  water vapour density from ERAINT 2-m temperature based on the ideal
  gas law (\code{Rho}).

  These indices were generated with the functions \code{\link{evaporation}},
  \code{\link{H2O}}, \code{\link{distill.index}},
  and \code{\link{moistflux}}.

  The time series objects \code{Es}, \code{H2O}, \code{Qs}, and
  \code{Rho} hold three different series: (is=1) the sum between
  50S-50N, (is=2) 90S-50S, and (is=3) 50N-90N.

  \code{Parea} contains nine series: Precipitation area,Fraction
  precipitation area, Precipitation area over ocean, Precipitation area
  over land, wet-day mean preciptation, total precipitation, total
  precipitation over oceans, total precipitation over land, and number
  of valid data. 
  
  \code{prtot.merra} and \code{prtot.eraint} are zoo objects that hold the mean total daily precipitation (unit=kilo-tons) for the latitude bands 50S-50N, 90S-50S, and 50N-90N, respectivey.
}
\usage{
data(Es)
data(H2O)
data(Parea)
data(Parea.merra)
data(Parea.merra.glob)
data(Parea.month)
data(Parea.merra.month)
data(Parea.eraint.month)
data(Parea.cmip5)
data(Qs)
data(Rho)
data(prtot.merra)
data(prtot.eraint)
data(raingaugedata.USA)
data(raingaugedata.World)
}