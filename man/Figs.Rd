\name{Fig1}
\alias{Fig1}
\alias{Fig2}
\alias{Fig3}
\alias{Fig4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Figures
%%  ~~function to do ... ~~
}
\description{Functions plotting the main figures in the paper
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Fig1(path = "~/TRMM", pattern = "3B42_Daily", param = "precipitation", nx = 1440, ny = 400, a = 6378, x0 = 1)
Fig2()
Fig3(col=c('black',rgb(0.5,0.5,0.5,0.5),'red',rgb(0.5,0,0,0.5),rgb(0,0,0.5,0.5)))
Fig4()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
%%     ~~Describe \code{path} here~~
}
  \item{pattern}{
%%     ~~Describe \code{pattern} here~~
}
  \item{param}{
%%     ~~Describe \code{varid} here~~
}
  \item{nx}{
%%     ~~Describe \code{nx} here~~
}
  \item{ny}{
%%     ~~Describe \code{ny} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{x0}{
%%     ~~Describe \code{x0} here~~
}
  \item{col}{colour scheme for the curves see \code{\link{plot}}.}
}

\author{Rasmus.Benestad-at-met.no
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
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
    nx = 1440, ny = 400, a = 6378, x0 = 0.1) 
{
    require(ncdf4)
    fname <- list.files(path = path, pattern = pattern, full.names = TRUE)[1]
    print(fname)
    ncid <- nc_open(fname)
    x <- ncvar_get(ncid, varid)
    x <- t(x)
    x[x < x0] <- NA
    x[x >= x0] <- 1
    dim(x) <- c(1, nx * ny)
    tim <- ncatt_get(ncid, varid = 0, attname = "EndDate")
    time <- tim$value
    lon <- ncvar_get(ncid, "lon")
    lat <- ncvar_get(ncid, "lat")
    nc_close(ncid)
    attr(x, "longitude") <- lon
    attr(x, "latitude") <- lat
    x <- as.field(zoo(x, order.by = time), lon = lon, lat = lat, 
        param = "precip", unit = "NA", src = "TRMM")
    map(x, type = "fill", colbar = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
