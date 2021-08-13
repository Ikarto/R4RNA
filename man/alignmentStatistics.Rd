\name{Alignment Statistics}
\alias{baseConservation}
\alias{basepairCanonical}
\alias{basepairConservation}
\alias{basepairCovariation}
\alias{helixCanonical}
\alias{helixConservation}
\alias{helixCovariation}
\alias{alignmentConservation}
\alias{alignmentCovariation}
\alias{alignmentCanonical}
\alias{alignmentPercentGaps}

\title{
	Compute statistics for a multiple sequence alignments
}
\description{
	Functions to compute covariation, percent identity conservation, and percent
	canonical basepairs given a multiple sequence alignment and optionally a
	secondary structure.  Statistics can be computed for a single base,
	basepair, helix or entire alignment.
}
\usage{
	baseConservation(msa, pos)

	basepairConservation(msa, pos.5p, pos.3p)
	basepairCovariation(msa, pos.5p, pos.3p) 
	basepairCanonical(msa, pos.5p, pos.3p)

	helixConservation(helix, msa)
	helixCovariation(helix, msa)
	helixCanonical(helix, msa) 

	alignmentConservation(msa)
	alignmentCovariation(msa, helix)
	alignmentCanonical(msa, helix)
	
	alignmentPercentGaps(msa)
}
\arguments{
	\item{helix}{
		A helix data.frame
	}
	\item{msa}{
		A multiple sequence alignment, such as those returned by
		\code{\link{readFasta}}
	}
	\item{pos, pos.5p, pos.3p}{
		Positions of bases or basepairs for which statistics shall be calculated
		for.
	}
}

\details{
	Conservation values have a range of [0, 1], where 0 is the absence of
	primary sequence conservation (all bases different), and 1 is full
	primary sequence conservation (all bases identical).
	
	Canonical values have a range of [0, 1], where 0 is a complete lack of
	basepair potential, and 1 indicates that all basepairs are valid

	Covariation values have a range of [-2, 2], where -2 is a complete lack of
	basepair potential and sequence conservation, 0 is complete sequence
	conservation regardless of basepairing potential, and 2 is a complete lack
	of sequence conservation but maintaining full basepair potential.

	\code{helix} values are average of base/basepair values, and the
	\code{alignment} values are averages of helices or all columns depending
	on whether the \code{helix} argument is required.	

	\code{alignmentPercentGaps} simply returns the percentage of nucleotides
	that are gaps in a sequence for each sequence of the alignment.
}

\value{
	\code{baseConservation}, \code{basepairConservation},
	\code{basepairCovariation}, \code{basepairCanonical},
	\code{alignmentConservation}, \code{alignmentCovariation}, and
	\code{alignmentCanonical}
	return a single decimal value.

	\code{helixConservation}, \code{helixCovariation}, \code{helixCanonical}
	return a list of values whose length equals the number of rows in \code{helix}.

	\code{alignmentPercentGaps} returns a list of values whose length equals
	the number of sequences in the multiple sequence alignment.
}

\author{
	Jeff Proctor, Daniel Lai
}

\examples{
	data(helix)
	
	baseConservation(fasta, 9)

	basepairConservation(fasta, 9, 18)
	basepairCovariation(fasta, 9, 18) 
	basepairCanonical(fasta, 9, 18)

	helixConservation(helix, fasta)
	helixCovariation(helix, fasta)
	helixCanonical(helix, fasta)

	alignmentConservation(fasta)
	alignmentCovariation(fasta, helix)
	alignmentCanonical(fasta, helix)

	alignmentPercentGaps(fasta)
}
\keyword{math}
