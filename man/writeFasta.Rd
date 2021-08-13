\name{Write FASTA}
\alias{writeFasta}
\title{
	Writes out a FASTA format file from a list of named characters
}
\description{
	Writes out a FASTA format file from a list of named characters, where the
	sequences are from the elements, and the descriptions are form the names.
	Does not attempt to break the sequence into multiple lines.
}
\usage{
	writeFasta(msa, file = stdout(), wrap = NA)
}
\arguments{
	\item{msa}{
		A list of characters representing each sequence, and names for each
		element containing the description of each sequence.  Defaults to the
		console.
	}
	\item{file}{
		A character string pointing to the path of a file, or a connection.
	}
	\item{wrap}{
		An integer to determine the number of characters in each row before
		the line wraps/breaks.  If NA, then no wrapping will occur.
	}
}

\author{
	Daniel Lai
}

\examples{
	fasta <- c(sequence = "AAAAACCCCCUUUUU", structure = "(((((.....)))))")
	writeFasta(fasta)
}
\keyword{IO}
\keyword{file}
