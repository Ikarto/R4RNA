\name{Read Structure File}
\alias{readHelix}
\alias{readConnect}
\alias{readVienna}
\alias{readBpseq}

\title{
	Read secondary structure file
}
\description{
	Reads in secondary structure text files into a helix data frame.
}
\usage{
	readHelix(file)
	readConnect(file)
	readVienna(file, palette = NA)
	readBpseq(file)
}
\arguments{
	\item{file}{
		A text file in connect format, see details for format specifications.
	}
	\item{palette}{
		Used to colour basepairs by bracket type.  See
		\code{\link{viennaToHelix}} for more details.
	}

}
\details{
	\bold{Helix}: Files start with a header line beginning with # followed by
	the sequence length, followed by a four-column tab-delimited table (with
	column names), where each row corresponds to a helix in the structure.  The
	four columns are \emph{i} and \emph{j} for the left-most and right-most
	basepair positions respectively, the \emph{length} of the helix (converging
	inwards from i and j, and finally an arbitrary \emph{value} assigned to the
	helix.
	
	\bold{Vienna}: Dot-bracket notation from Vienna package programs, where each
	structure consists of matched brackets for basepairs and periods for unbased
	pairs.  Valid brackets are (, {, [, <, A, B, C, D matched with ), }, ], >,
	a, b, c, d, respectively.  An energy value can be appended to
	the end of any dot-bracket structure.  The function will accept slight
	variations of the format, including those with FASTA-like headers (in which
	case line breaks are allows), and those without FASTA-like headers (in which
	case line breaks are NOT allowed), with both types allowing for a preceding
	(NOT following) nucleotide sequence for the structure.  Multiple entries
	\emph{of the same length} may be in a single file, which will be returned
	as a single helix structure, with respectively energy values (if specified).
	
	\bold{Connect}: Output from mfold and other programs, this format is
	expected to be a text file beginning with a header line that starts with the
	sequence length, with an optional Energy/dG value, followed by a six-column
	tab-delimited table where columns 1 and 5 denote the position that are
	basepaired (unpaired when column 5 is 0).  Other columns are ignored, but
	for completeness, column 2 is the nucleotide, column 3 and 4 are the
	positions of the bases left and right of the base specified in column 1
	respectively (with 0 denoting non-existance), and column 6 a copy of column
	1. Multiple entries \emph{of the same length} may be in a single file, which
	will be returned as a single helix structure. All helices will be assigned
	the energy value extracted from their respective structure header lines.
	
	\bold{Bpseq}: Format used by the Gutell Lab's Comparative RNA Website.  The
	file may optionally begin several header lines (e.g. Filename,
	Organism, Accession, etc.), followed by a 3-column tab-delimited table for
	the structure, where column 1 is the base position, base 2 is the nucleotide
	base, and column 3 is the paired position (0 if unpaired).  Certain pieces
	of header information will be parsed and returned as attributes of the
	output data frame.  Multiple structures can be within a single file,
	returned as a single helix data frame, with attributes set to those of the
	first entry.

}
\value{
	Returns a helix format data frame.
}
\author{
	Daniel Lai, Jeff Proctor
}

\examples{
	file <- system.file("extdata", "helix.txt", package = "R4RNA")
	helix <- readHelix(file)
	head(helix)

	file <- system.file("extdata", "connect.txt", package = "R4RNA")
	connect <- readConnect(file)
	head(connect)
	message("Note connect data assigns structure energy level to all basepairs")

	file <- system.file("extdata", "vienna.txt", package = "R4RNA")
	vienna <- readVienna(file)
	head(vienna)
	message("Note vienna data assigns structure energy level to all basepairs")

	file <- system.file("extdata", "bpseq.txt", package = "R4RNA")
	bpseq <- readBpseq(file)
	head(bpseq)
	message("Note bpseq data has no value assigned to basepairs")
}
\keyword{IO}
\keyword{file}
