\name{timeseries}
\alias{timeseries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Time series
%%  ~~function to do ... ~~
}
\description{Plots the time series 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
timeseries(x, is = 1, denominator='total area', ylab=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a time series object
%%     ~~Describe \code{x} here~~
}
  \item{is}{index - which series to plot
%%     ~~Describe \code{is} here~~
}
\item{denominator}{Denominator}
\item{ylab}{see \code{link{plot}}}
  \item{\dots}{arguments passed on to \code{\link{plot}}
%%     ~~Describe \code{\dots} here~~
}
}

\author{Rasmus.Benestad@met.no
%%  ~~who you are~~
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
function (x, is = 2, ...) 
{
    par(bty = "n")
    y <- x[, is]
    plot(y, main = attr(x, "longname")[is], plot.type = "single", 
        col = c("red", "blue", "grey"), sub = paste("TMRR ", 
            paste(abs(range(lat(x))), collapse = "S-"), "N: ", 
            min(year(y)), " - ", max(year(y)), sep = ""), ylab = attr(x, 
            "unit")[is], xlab = "", ylim = c(0, max(y)), ...)
    for (i in 1:length(is)) {
        lines(trend(y[, i]))
        t1 <- trend(y[, i])[1]
        t2 <- trend(y[, i])[length(index(y))]
        lines(c(as.Date("1990-01-01"), index(y)[1]), rep(t1, 
            2), col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 2)
        lines(c(as.Date("1990-01-01"), index(y)[length(index(y))]), 
            rep(t2, 2), col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 2)
        text(as.Date("1997-08-01"), t1, round(t1), cex = 0.7)
        text(as.Date("1997-08-01"), t2, round(t2), cex = 0.7)
    }
    grid()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ timeseries plot }
\keyword{ trend analyss }% __ONLY ONE__ keyword per line
