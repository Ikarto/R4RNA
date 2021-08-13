\name{Read FASTA}
\alias{readFasta}
\title{
	Read FASTA format multiple sequence alignment files.
}
\description{
	Reads in FASTA format multiple sequence text files into a list of
	named characters, with names derived from the description line of each
	FASTA entry.
}
\usage{
	readFasta(file, filter = FALSE)
}
\arguments{
	\item{file}{
		FASTA format file containing at least one sequence, where each sequence
		has a description line beginning with the > character.
	}
	\item{filter}{
		When true, filters out any sequences with any other characters besides:
		A, C, G, T, U, and - for gaps.  Also converts all letter to uppercase
		and all T's to U's, and N's to -'s.
	}
}
\value{
	Returns an array of named characters, each element a sequence read, with
	the description as its name.
}
\author{
	Daniel Lai
}

\examples{
	file <- system.file("extdata", "fasta.txt", package = "R4RNA")
	fasta <- readFasta(file)
	head(fasta)
}
\keyword{IO}
\keyword{file}
