\name{Convert Helix Formats}
\alias{viennaToHelix}
\alias{helixToVienna}
\alias{helixToConnect}
\alias{helixToBpseq}
\title{
	Convert helix structures to and from other formats
}
\description{
	Converts dot bracket vienna format to and from helix format.  It should be
	noted that the allows structures of vienna is a subset of those allowed
	in the helix format.  Thus, conversion from vienna to helix will yield the
	identical structure, while conversion from helix to vienna may result in the
	loss of certain basepairs (mainly those that are conflicting).  Pseudoknots
	are supported in both directions of conversion with limitations.
}
\usage{
	viennaToHelix(vienna, value = NA, palette = NA)
	helixToVienna(helix)
	helixToConnect(helix)
	helixToBpseq(helix)
}
\arguments{
	\item{vienna}{
		A string containing \emph{only} a vienna dot bracket structure, with
		balanced brackets.  Allowable brackets are (, <, [, \{, A, B, C, and D
		(where upper-case alphabets are paired with lower-case alphabets).
	}
	\item{value}{
		A numerical value to assign to all helices.
	}
	\item{palette}{
		A list of colour names for up to 8 colours that will be used to colour
		brackets of type (, <, [, \{, A, B, C, and D, respectively.
	}
	\item{helix}{
		A helix data.frame.
	}
}
\details{
	\code{viennaToHelix} will ignore any non dot-bracket characters prior
	to parsing, so the resultant length will be shorter than expected if invalid
	characters are included.

	If the colour palette is less than the number of supported brackets, it
	will simply cycle through the list.  To explicitly prevent the colouring/
	display of specific bracket type, colour it \dQuote{NA}.

	For \code{helixToVienna}, pseudoknotted basepairs will be assigned
	different bracket types.  As there are only 8 supported bracket types,
	any basepair pseudonotted deeper than 8 levels will be excluded from the
	output.  Additionally, vienna format is unable to respresent conflicting
	basepairs, so conflicting basepairs will also be excluded.  For both types
	of exclusion, those at the bottom of the helix data.frame will always be
	excluded in favour of keeping helices higher on the data.frame table.	

	\code{helixToConnect} and \code{helixToBpseq} will convert a
	\emph{non-conflicting} helix data.frame into connect or bpseq format
	repsectively, provided the helix structure has a \dQuote{sequence} attribute
	containing a single nucleotide sequence of the structure.
}
\author{
	Daniel Lai
}
\examples{
	# viennaToHelix demonstrating ALL valid bracket symbols
	dot_bracket <- ".....(<[{.....ABCD.....}]>).....dcba....."
	parsed <- viennaToHelix(dot_bracket, -31.5)
	print(parsed)

	vienna <- helixToVienna(parsed)
	print(vienna)

	# Colouring the brackets by bracket type
	colour <- c("red", "orange", "yellow", "green", "lightblue", "blue", "purple", "black")
	double.rainbow <- viennaToHelix(dot_bracket, 0, colour)
	plotHelix(double.rainbow)
}
\keyword{manip}
