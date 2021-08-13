\name{Coerce to Helix}
\alias{as.helix}
\alias{is.helix}
\alias{parseBracket}
\title{
	Coerce to a Helix Data Frame
}
\description{
	Functions to coerce a structure into a helix data frame, and to check
	whether a structure is a valid helix data frame.  A helix data frame is
	a data frame, so any structure coercible into a data.frame can become a
	helix data frame.
}
\usage{
	as.helix(x, length)
	is.helix(x)
}
\arguments{
	\item{x}{
		Structure to coerce.  Should be a structure coercible into a standard
		R data.frame structure for \code{\link{as.helix}}. Should be a string for
		\code{\link{parseBracket}}.  May be anything for \code{\link{is.helix}}.
	}
	\item{length}{
		The length of the RNA sequence containing the helices.
	}
}
\details{
	\code{\link{as.helix}} takes in a data.frame and coerces it into a 
	helix data frame acceptable by other R4RNA functions.  This mainly involves
	setting specific column names and casting to specific types.
}
\author{
	Daniel Lai
}
\examples{
	# Not a valid helix data frame
	helix <- data.frame(c(1, 2, 3), seq(10, 20, length.out = 3), 5, runif(3))
	is.helix(helix)
	warnings()
	
	# Formatted into a helix data frame
	helix <- as.helix(helix)
	is.helix(helix)
}
\keyword{manip}
