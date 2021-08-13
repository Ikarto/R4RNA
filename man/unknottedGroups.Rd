\name{Find Unknotted Groups}
\alias{unknottedGroups}
\title{
	Partition basepairs into unknotted groups
}
\description{
	Breaks down input helices into basepairs, and assigns each basepair to a
	numbered group such that basepairs in each group are non-pseudoknotted
	relative to all other basepairs within the same group.

	The algorithm is greedy and thus will \emph{not} find the best combination
	of basepairs to minimize the number of groups.
}
\usage{
	unknottedGroups(helix) 
}
\arguments{
	\item{helix}{
		A helix data.frame.
	}
}
\value{
	An array of integers dictating the groups of each helix.  Will only
	correspond to the input helix structure if the input had helices of length 1
	(e.g. output of \code{\link{expandHelix}}).
}
\seealso{
	\code{\link{colourByUnknottedGroups}}
	\code{\link{expandHelix}}
}
\author{
	Daniel Lai
}
\examples{
	data(helix)
	known$group <- unknottedGroups(known)
	print(known)
}
\keyword{aplot}
