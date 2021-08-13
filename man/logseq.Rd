\name{Log10 Space Operations}
\alias{logseq}
\alias{logfloor}
\alias{logceiling}

\title{
	Log base 10 sequence, floor and ceiling
}
\description{
	Sequence, floor and ceiling operations in log 10 space.
}
\usage{
	logseq(from, to, length.out)

	logfloor(x)

	logceiling(x)
}
\arguments{
	\item{from, to}{
		Positive non-zero values to start and end sequence, respectively.
	}
	\item{length.out}{
		The number of elements the resulting sequence should containg. If
		absent, function will attempt to generate numbers factors of 10 apart.
	}
	\item{x}{
		A value to round.
	}
}
\value{
	\code{logseq} returns an array numbers evenly distanced in log10-space.
	
	\code{logfloor} and \code{logceiling} return a value that is 10 raised
	to an integer number.
}

\author{
	Daniel Lai
}

\examples{
	logseq(1e-10, 1e3)
	logseq(1e-10, 1e3, length.out = 10)
	logceiling(2.13e-6)
	logfloor(2.13e-6)
}
\keyword{color}
