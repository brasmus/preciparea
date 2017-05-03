\name{evaporation}
\alias{evaporation}
\alias{moistflux}
\alias{columnH2O}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaporation estimation
%%  ~~function to do ... ~~
}
\description{Processing routines applied to ERAINT reanalysis to estiamte global evaporation and atmospheric moisture content
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
evaporation(fname = "~/Downloads/ERAINT_t2m.nc", R.s = 18.02)
moistflux(fname = "~/Downloads/ERAINT_Qs.nc")
columnH2O(fname = "~/Downloads/ERAINT_Q-columntotal.nc")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{File name
%%     ~~Describe \code{fname} here~~
}
  \item{R.s}{specific gas constant for water vapour
%%     ~~Describe \code{R.s} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{time series objects \code{\link{Es}}, \code{\link{Rho}}, \code{\link{Qs}}, and \code{\link{H2O}}. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
function (fname = "~/Downloads/ERAINT_t2m.nc", R.s = 18.02) 
{
    t2m <- retrieve(fname)
    es <- C.C.eq(mask(t2m, land = TRUE))
    rho <- es/(R.s * (t2m + 273.15))
    es.50S50N <- aggregate.area(subset(es, is = list(lon = range(lon(es)), 
        lat = c(-50, 50))), FUN = "mean")
    es.90S50S <- aggregate.area(subset(es, is = list(lon = range(lon(es)), 
        lat = c(-90, -50))), FUN = "mean")
    es.50N90N <- aggregate.area(subset(es, is = list(lon = range(lon(es)), 
        lat = c(50, 90))), FUN = "mean")
    rho.50S50N <- aggregate.area(subset(rho, is = list(lon = range(lon(es)), 
        lat = c(-50, 50))))
    rho.90S50S <- aggregate.area(subset(rho, is = list(lon = range(lon(es)), 
        lat = c(-90, -50))))
    rho.50N90N <- aggregate.area(subset(rho, is = list(lon = range(lon(es)), 
        lat = c(50, 90))))
    Es <- combine.stations(es.50S50N, es.90S50S, es.50N90N)
    Rho <- combine.stations(rho.50S50N, rho.90S50S, rho.50N90N)
    save(Es, file = "Es.rda")
    save(Rho, file = "Rho.rda")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ moisture }
\keyword{ evaporation }% __ONLY ONE__ keyword per line
