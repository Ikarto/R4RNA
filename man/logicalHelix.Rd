\name{Helix Type Filters}
\alias{isConflictingHelix}
\alias{isDuplicatingHelix}
\alias{isOverlappingHelix}
\title{
	Logical filters of helix by type
}
\description{
	Given a helix data frame, checks if helices are conflicting, duplicating,
	or overlapping, and returns an array of logicals.  See details for exact
	definition of the three types of events.
}
\usage{
	isConflictingHelix(helix, any = TRUE)
	isDuplicatingHelix(helix, any = TRUE)
	isOverlappingHelix(helix, query, any = FALSE)
}
\arguments{
	\item{helix}{
		A helix data frame
	}
	\item{any}{
		Relevant only for helices of length greater than 1.  If TRUE, then if
		any single basepair of the helix satisfies the condition, then the
		entire helix satisfies the condition.  If FALSE, then a helix only
		satisfies when all of its basepairs satisfy the condition.
	}
	\item{query}{
		For \code{isOverlappingHelix}, a helix data structure against which
		\code{helix} will be checked for overlap against.
	}
}
\details{
	Helices of length greater than 1 are internall expanded into basepairs of
	length 1, after which the following conditions are evaluated:

	A \bold{conflicting} basepair is one where at least one of its two positions
	is used by either end of another basepair.

	A \bold{duplicating} basepair is one where both of its positions are used
	by both ends of another basepair.

	An \bold{overlapping} basepair is one in \code{helix} where both of its
	positions are used by both ends of another basepair in the \code{query}
	structure.

	In the case of \emph{conflicting} and \emph{duplicating} basepairs, for a
	set of basepairs that satisfies this condition, the basepair situation
	highest on the data frame will be exempt from the condition. i.e.  Say 5
	basepairs are all duplicates of each other, the top 1 will return FALSE,
	while the bottom 4 will return TRUE.  This assumes some significant meaning
	to the ordering of rows prior to using this function.  This is to be used
	with \code{which} to filter out basepairs that satisfy these conditions,
	leaving a set of basepairs free of these events.

	If the original input had helices greater than length 1, then after applying
	all of the above, a conditional is assigned to each of the original input
	helices according to the rule set by \code{any}.
}
\value{
	Returns an array of logicals corresponding to each row of \code{helix},
	giving the conditional status of the helix.
}
\author{
	Daniel Lai
}

\examples{
	data(helix)

	conflicting <- isConflictingHelix(helix)
	duplicating <- isDuplicatingHelix(helix)

	# Nonsensical covariation plot
	plotCovariance(fasta, helix)

	# Plot nonconflicting helices
	plotCovariance(fasta, helix[(!conflicting & !duplicating), ])

	# Similar result
	plotCovariance(fasta, helix, conflict.filter = "lightgrey")
}
\keyword{logic}

