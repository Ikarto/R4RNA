\name{Colour Helices}
\alias{colourByCount}
\alias{colourByValue}
\alias{colourByConservation}
\alias{colourByCovariation}
\alias{colourByCanonical}
\alias{colourByBasepairFrequency}
\alias{colourByUnknottedGroups}
\alias{defaultPalette}

\title{
	Assign colours to helices
}
\description{
	Functions to generate colours for helices by various rules, including
	integer counts, value ranges, percent identity covariation, conservation,
	percentage canonical basepair, basepair frequency, and non-pseudoknotted
	groups.
}
\usage{
	colourByCount(helix, cols, counts, get = FALSE)
	colourByValue(helix, cols, breaks, get = FALSE,
		log = FALSE, include.lowest = TRUE, ...)
	colourByBasepairFrequency(helix, cols, get = TRUE)
	colourByUnknottedGroups(helix, cols, get = TRUE)
	colourByCovariation(helix, msa, cols, get = FALSE)
	colourByConservation(helix, msa, cols, get = FALSE)
	colourByCanonical(helix, msa, cols, get = FALSE)
	defaultPalette()
}
\arguments{
	\item{helix}{
		A helix data frame to be coloured.
	}
	\item{cols}{
		An array of characters (or numbers) representing a set of colours to
		colour \code{helix} with.  When missing, a default set of colours from
		defaultPalette() will be used.  Valid input include hex codes, colour
		names from the \code{colours} function, and integer numbers.  The
		colours will be interpreted as being from best to worst.
	}
	\item{counts}{
		An array of integers the same length as \code{cols}, dictating the
		number of times each corresponding colour should be used.  When missing,
		the function will divide the number of helices evenly over each of the
		colours available.
	}
	\item{breaks}{
		An integer number of intervals to break the \sQuote{value} column of
		\code{helix} into, or a list of numbers defining the interval breaks.
		If missing, the range of \sQuote{helix$value} will automatically be
		split evenly into intervals for each colour available.
	}
	\item{get}{
		If TRUE, returns the input \code{helix} with a \code{col} column, else
		simply returns an array of colours the same length as the number of row
		in \code{helix}.  The exceptions are \code{colourByBasepairFrequency}
		and \code{colourByUnknottedGroups} which will return a different helix
		if TRUE, and a list of colours that will not match the input helix if
		FALSE.
	}
	\item{log}{
		If TRUE, will breaks values into even log10 space intervals, useful
		when values are p-values.
	}
	\item{include.lowest}{
		Whether the lowest interval should include the lowest value, passed to
		\code{\link{cut}}
	}
	\item{...}{
		Additional arguments passed to \code{\link{cut}}, potentially useful
		ones include \code{right} (whether intervals should be inclusive on the
		right or left) and \code{dig.lab} (number of digits in interval labels).
	}
	\item{msa}{
		A multiple sequence alignment, such as those returned by
		\code{\link{readFasta}}
	}
}
\details{
	\code{colourByCount} assigns colours indepenent of the helix input's value
	column, and instead operates over the number of helices (i.e. rows).
	
	\code{colourByValue} uses \code{\link{cut}} to assign each of the helices
	to an interval based on its value.
	
	\code{colourByCovariation}, \code{colourByConservation} , and
	\code{colourByCanonical}, colour helices according to compensatory mutations
	(or covariation), percentage identity conservation, and percentage canonical
	basepair repsectively, relative to the multiple sequence alignment provided.

	\code{colourByBasepairFrequency} colours each basepair according to the
	number of times it appear in the input, regardless of its value.

	\code{colourByUnknottedGroups} greedily partitions the basepairs into non-
	pseudoknotted groups, and assigns a colour to each.
}

\value{
	All \dQuote{colourBy} functions return a list of colours when \code{get =
	FALSE}, and a helix with a \code{col} column if \code{get = TRUE}.  In both
	bases, the returned object has attributes \dQuote{legend} and \dQuote{fill},
	showing the mapping between interval (in legend) and colour (in fill),
	which can as eponymous arguments \code{\link{legend}}.
	
	\code{defaultPalette} returns the default list of colours.
}
\seealso{
	\code{\link{plotHelix}}
	\code{\link{logseq}}
	\code{\link{basepairFrequency}}
	\code{\link{unknottedGroups}}
}

\author{
	Daniel Lai
}

\examples{
	data(helix)

	known$col <- colourByCount(known)
	plotHelix(known)

	plotHelix(colourByValue(helix, log = TRUE, get = TRUE))

	cov <- colourByCovariation(known, fasta, get = TRUE)
	plotCovariance(fasta, cov)
	legend("topleft", legend = attr(cov, "legend"),
		fill = attr(cov, "fill"), title = "Covariation")

}
\keyword{color}
