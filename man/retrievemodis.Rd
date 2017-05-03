\name{retrievemodis}
\alias{retrievemodis}
\alias{retrieveMERRA}
\alias{TRMM2monthly}
\alias{monthP2area}
\alias{PareaCMIP}
\alias{globalsum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Retrieve MODIS
%%  ~~function to do ... ~~
}
\description{Function that reads all the TRMM or MERRA data and distills indices.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
retrievemodis(path = "~/TRMM", pattern = "3B42_Daily", varid = "precipitation", nx = 1440, ny = 400, a = 6378, x0 = 1)
retrieveMERRA(path='~/MERRA',
              pattern='MERRA2_',reg50s50n=TRUE,
              param=c('PRECCUCORR','PRECLSCORR','PRECSNOCORR'),
              x0 = 1) 
TRMM2monthly(path='~/TRMM',
             pattern='3B42_Daily',
             param='precipitation',yr1=1998,yr2=2016,
             nx=1440,ny=400)          
monthP2area(fname='~/TRMM-monthly/trmm-precip-mon.nc',param='precip',x0=1,dy=0.1)
PareaCMIP(path='CMIP5.monthly/rcp45/',pattern='pr_Amon_ens_',param='pr',x0=4,dy=0.1)
globalsum(fname,FUN='sum',dlat=10,dy=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{path where the TRMM data is kept
%%     ~~Describe \code{path} here~~
}
  \item{pattern}{Text pattern of the TRMM data
%%     ~~Describe \code{pattern} here~~
}
  \item{param}{variable name
%%     ~~Describe \code{varid} here~~
}
  \item{nx}{number of grid boxes along the zonal direction
%%     ~~Describe \code{nx} here~~
}
  \item{ny}{number of grid boxes along the meridional direction
%%     ~~Describe \code{ny} here~~
}
  \item{a}{radius of the erth (km)
%%     ~~Describe \code{a} here~~
}
  \item{x0}{threshold to distinguish between dry and wet grid box.
%%     ~~Describe \code{x0} here~~
}
\item{fname}{Name of the netCDF file with monthly mean precipitation}
\item{dlat}{the increments in latitudes}
\item{dy}{smaller than the increments in latitudes}
\item{reg50s50n=TRUE}{FALSE for global; TRUE for 50S--50N}
}

\author{Rasmus.Benestad@met.no
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
## Retrieving the data:
## https://disc2.gesdisc.eosdis.nasa.gov/data/TRMM_L3/TRMM_3B42_Daily.7/
## https://disc.sci.gsfc.nasa.gov/recipes/?q=recipes/How-to-Download-Data-Files-from-HTTP-Service-with-wget
## https://disc2.gesdisc.eosdis.nasa.gov/data/TRMM_L3/TRMM_3B42_Daily.7/doc/TRMM_Readme_v3.pdf

## Script for downloading thr TRMM data
#!/bin/bash

for year in {1998..2016}
do
  for month in {1..12}
  do
  wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies 
       -r -c -nH -nd -np -A nc4 "https://disc2.gesdisc.eosdis.nasa.gov/data/TRMM_L3/TRMM_3B42_Daily.7/$year/$month"
  done
done

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (path = "~/TRMM", pattern = "3B42_Daily", varid = "precipitation", 
    nx = 1440, ny = 400, a = 6378, x0 = 1e-10) 
{
    require(ncdf4)
    fnames <- list.files(path = path, pattern = pattern, full.names = TRUE)
    print(fnames)
    nt <- length(fnames)
    Ap <- rep(NA, nt)
    Af <- Ap
    Am <- Ap
    Pt <- Ap
    mu <- Ap
    Pm <- Ap
    nv <- Ap
    Al <- Ap
    Pl <- Ap
    t <- rep(NA, nt)
    for (i in 1:nt) {
        fname <- fnames[i]
        print(fname)
        ncid <- nc_open(fname)
        x <- ncvar_get(ncid, varid)
        x <- t(x)
        dim(x) <- c(1, nx * ny)
        tim <- ncatt_get(ncid, varid = 0, attname = "EndDate")
        t[i] <- tim$value
        lon <- ncvar_get(ncid, "lon")
        lat <- ncvar_get(ncid, "lat")
        nc_close(ncid)
        attr(x, "longitude") <- lon
        attr(x, "latitude") <- lat
        z <- x
        zc <- x
        if (i == 1) {
            Lon <- rep(lon, ny)
            Lat <- sort(rep(lat, nx))
            aweights <- cos(pi * Lat/180)
            aweights <- aweights/mean(aweights)
            A.tot <- sum(aweights)
            test <- subset(retrieve("air.mon.mean.nc"), it = 1)
            test <- regrid(test, is = list(lon = lon, lat = lat))
            ocean <- mask(test, land = TRUE)
            land <- mask(test, land = FALSE)
            ocean <- coredata(ocean)
            ocean[is.finite(ocean)] <- 0
            land <- coredata(land)
            land[is.finite(land)] <- 0
        }
        zl <- x + land
        zo <- x + ocean
        Pm[i] <- sum(zo * aweights, na.rm = TRUE)
        zo[zo > x0] <- 1
        zo[zo <= x0] <- 0
        Am[i] <- sum(zo * aweights, na.rm = TRUE)
        Pl[i] <- sum(zl * aweights, na.rm = TRUE)
        zl[zl > x0] <- 1
        zl[zl <= x0] <- 0
        Al[i] <- sum(zl * aweights, na.rm = TRUE)
        zc[zc > x0] <- 1
        zc[zc <= x0] <- 0
        Ap[i] <- sum(zc * aweights, na.rm = TRUE)
        Pt[i] <- sum(x * aweights, na.rm = TRUE)
        z[z <- x0] <- NA
        mu[i] <- mean(z * aweights, na.rm = TRUE)
        nv[i] <- sum(is.finite(x))
    }
    Fp <- round(100 * Ap/A.tot, 2)
    X <- cbind(a * Ap, Fp, a * Am, a * Al, mu, Pt, Pm, Pl, nv)
    colnames(X) <- c("A.precip", "A.precip.frac", "A.precip.marine", 
        "A.precip.land", "\n                   mu", "P.tot", 
        "P.marine", "P.land", "nv")
    Y <- zoo(X, order.by = as.Date(t))
    attr(Y, "longitudes") <- range(lon)
    attr(Y, "latitudes") <- range(lat)
    attr(Y, "unit") <- c("km^2", "\%", "km^2", "km^2", "mm", "mm", 
        "mm", "mm", "count")
    attr(Y, "longname") <- c("Precipitation area", "Fraction precipitation area", 
        "Precipitation area over ocean", "Precipitation area over land", 
        "wet-day mean preciptation", "total precipitation", "total precipitation over oceans", 
        "total precipitation over land", "number of valid data")
    attr(Y, "src") <- "TRMM"
    attr(Y, "precip.threshold") <- x0
    attr(Y, "total area") <- A.tot
    attr(Y, "url") = "http://disc2.gesdisc.eosdis.nasa.gov/data/TRMM_L3/TRMM_3B42_Daily.7/"
    attr(Y, "history") <- "trmm-precip.R"
    attr(Y, "author") <- "rasmus.benestad@met.no"
    return(Y)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ TRMM }
\keyword{ precipitation }% __ONLY ONE__ keyword per line
