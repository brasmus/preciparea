\name{PDF}
\alias{PDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Probability density function and histogram
%%  ~~function to do ... ~~
}
\description{A graph showing the marginal distribution of the time series
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PDF(x, is = 1, plot=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a time series
%%     ~~Describe \code{x} here~~
}
  \item{is}{index for which series to show
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

\seealso{\code{\link{AC}}, \code{\link{timeseries}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, is = 1, ...) 
{
    z <- coredata(trend(x[, is], result = "residual"))
    hist(z, col = "grey", freq = FALSE, main = "", ylab = "", 
        xlab = paste(attr(x, "longname")[is], " (", attr(x, "unit")[is], 
            ")", sep = ""), ...)
    x <- seq(-max(abs(z)), max(abs(x)), length = 100)
    lines(x, dnorm(x, mean = mean(z), sd = sd(z)), lwd = 3, col = rgb(1, 
        0, 0, 0.4))
    grid()
    par(new = TRUE, fig = c(0.75, 0.98, 0.75, 0.98), mar = rep(0, 
        4), cex.axis = 0.7)
    qqnorm(z, main = "")
    qqline(z, col = "red")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ histogram }
\keyword{ statistical distribution }% __ONLY ONE__ keyword per line
