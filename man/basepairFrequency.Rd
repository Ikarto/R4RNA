\name{Basepair Frequency}
\alias{basepairFrequency}
\title{
	Calculates the frequency of each basepair
}
\description{
	Calculates the frequency of each basepair in a given helix structure.
	Internally, breaks helices into basepairs, and returns a structure of
	unique basepairs, where the values is its frequency, regardless of original
	value.
}
\usage{
	basepairFrequency(helix) 
}
\arguments{
	\item{helix}{
		A helix data.frame
	}
}
\value{
	A helix data.frame of unique basepairs of length 1, with the frequency
	of appearance as its value, sorted by decreasing value.
}
\seealso{
	\code{\link{colourByBasepairFrequency}}
}
\author{
	Daniel Lai
}
\examples{
	data(helix)
	basepairFrequency(helix)
}
\keyword{aplot}
