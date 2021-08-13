\name{Structure Mismatch Score}
\alias{structureMismatchScore}
\title{
	Scores how a basepair structure fits a sequence
}
\description{
	Calculates a score that indicates how badly a set of basepairs (i.e. a
	secondary structure) fits with a sequence.  A perfect fit is a structure
	where all basepairs form valid basepairs (A:U, G:C, G:U, and equivalents)
	and has a score of 0.  Each basepair that forms a non-canonical pairing
	or pairs to gaps increases the score by 1, and each base-pair with a
	single-sided gap increases the score by 2.
}
\usage{
	structureMismatchScore(msa, helix)
}
\arguments{
	\item{msa}{
		An array of strings representing sequences of interest, typically the
		output from \code{\link{readFasta}}
	}
	\item{helix}{
		A helix data.frame
	}
}
\value{
	Returns an array of mismatch scores.
}
\author{
	Jeff Proctor, Daniel Lai
}
\examples{
	data(helix)
	mismatch <- structureMismatchScore(fasta, known)
	
	# Sort by increasing mismatch
	sorted_fasta <- fasta[order(mismatch)]
}
\keyword{math}
