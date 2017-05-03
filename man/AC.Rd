\name{AC}
\alias{AC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Annual cycle
%%  ~~function to do ... ~~
}
\description{Plots the mean annual cycle making use of the information embedded in this systematic and foced variation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AC(x, is = 1, plot=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a time series object
%%     ~~Describe \code{x} here~~
}
  \item{is}{index - which series to plot
%%     ~~Describe \code{is} here~~
}
\item{plot}{true: plot the results}
  \item{\dots}{arguments passed on to \code{\link{plot}}
%%     ~~Describe \code{\dots} here~~
}
}
\author{Rasmus.Benestad@met.no
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{PDF}}, \code{\link{timeseries}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, is = 2, ...) 
{
    y <- as.monthly(as.station(x[, is]), FUN = "mean")
    plot(aggregate(y, month, FUN = "mean"), col = c("red", "blue", 
        "grey"), lwd = 3, errorbar = FALSE, main = attr(x, "longname")[is], 
        ylab = attr(x, "unit")[is], xlab = "month", new = FALSE, 
        ...)
    grid()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ annual cycle }
\keyword{ systematic variations }% __ONLY ONE__ keyword per line
