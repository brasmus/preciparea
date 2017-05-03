\name{readGDCN}
\alias{readGDCN}
\alias{raingauges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rain gauges
%%  ~~function to do ... ~~
}
\description{Functions to read rain gauge data from GCDN and aggregate statistics.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
readGDCN(name, type = "PRCP")
raingauges(x0=1,tim=seq(as.Date('1960-01-01'),as.Date('2016-01-01'),by=1),
                       max.stations=NULL,
                       N.min=40*360,path="~/GDCN/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{file name
%%     ~~Describe \code{name} here~~
}
  \item{type}{variable
%%     ~~Describe \code{type} here~~
}
\item{x0}{threshold defining precipitation (mm/day)}
\item{tim}{time period}
\item{max.stations}{for testing with a smaller number of stations}
\item{N.min}{mnimum number of valid data points in the rain gauge records for inclusion}
\item{path}{path to the data files}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
function (name, type = "PRCP") 
{
    cnames <- c(paste("day.", 1:31, sep = ""), paste("flags.day.", 
        1:31, sep = ""))
    reshufle <- rep(0, 62)
    ii <- 0
    for (i in seq(1, 62, by = 2)) {
        ii <- ii + 1
        reshufle[i:(i + 1)] <- seq(ii, 62, by = 31)
    }
    cnames <- cnames[reshufle]
    x <- read.fwf(name, widths = c(3, 8, 4, 2, 4, rep(c(5, 2), 
        31)), col.names = c("country", "stnr", "year", "month", 
        "type", cnames))
    ipick <- is.element(x$type, type)
    if (sum(ipick) > 0) {
        dat <- as.matrix(x[ipick, seq(6, 67, by = 2)]) * 0.1
        dat[dat < -99] <- NA
        attr(dat, "Data source") <- "GDCN"
        attr(dat, "year") <- x$year[ipick]
        attr(dat, "month") <- x$month[ipick]
        attr(dat, "URL") <- "http://www.ncdc.noaa.gov/oa/climate/research/gdcn/gdcn.html"
        attr(dat, "Station_number") <- x$stnr[1]
        attr(dat, "Country_code") <- x$country[1]
        attr(dat, "original data name") <- name
        attr(dat, "flags") <- x[ipick, seq(6, 66, by = 2)]
        attr(dat, "history") <- "read with readGDCN - R-script."
        attr(dat, "Observations") <- rownames(table(x$type))
    }
    else dat <- NA
    invisible(dat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
