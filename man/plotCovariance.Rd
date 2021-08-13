\name{Covariation Plots}
\alias{plotCovariance}
\alias{plotDoubleCovariance}
\alias{plotOverlapCovariance}

\title{
	Plot nucleotide sequence coloured by covariance
}
\description{
	Given a multiple sequence alignment and a corresponding secondary structure,
	nucleotides in the sequence alignment will be coloured according to the
	basepairing and conservation status, where green is the most commonly
	observed valid basepair in the column, dark blue being valid covariation
	(i.e. mutation into another valid basepair), cyan is one-sided mutation that
	retains the basepair, and red is a mutation where the basepair has been
	lost.
}
\usage{
	plotCovariance(msa, helix, arcs = TRUE, add = FALSE, grid = FALSE,
		text = FALSE, legend = TRUE, species = 0, base.colour = FALSE, palette =
		NA, flip = FALSE, grid.col = "white", grid.lwd = 0, text.cex = 0.5,
		text.col = "white", text.font = 2, text.family = "sans", species.cex =
		0.5, species.col = "black", species.font = 2, species.family = "mono",
		shape = "circle", y = 0, conflict.lty = 2, conflict.col = NA,
		pad = c(0, 0, 0, 0), ...)
	plotDoubleCovariance(top.helix, bot.helix, top.msa, bot.msa = top.msa,
		add = FALSE, grid = FALSE, species = 0, legend = TRUE,
		pad = c(0, 0, 0, 0), ...)
	plotOverlapCovariance(predict.helix, known.helix, msa, bot.msa = TRUE,
		any = FALSE, miss = "black", add = FALSE, grid = FALSE, species = 0,
		legend = TRUE, pad = c(0, 0, 0, 0), ...)
}
\arguments{
	\item{msa, top.msa, bot.msa}{
		Multiple sequence alignment as an array of named characters, all of
		equal length.  Typically output of \code{\link{readFasta}}

		\code{top.msa} and \code{bot.msa} are specific to \code{top.helix} and
		\code{bot.helix} respectively, and may be set to NA to have no
		multiple sequence alignment at all.
	}
	\item{helix, top.helix, bot.helix, predict.helix, known.helix}{
		A helix data.frame with a structure corresponding to \code{msa},

		See \code{\link{plotDoubleHelix}} and \code{\link{plotOverlapHelix}} for
		detailed explanations of \code{top.helix}, \code{bot.helix},
		\code{predict.helix}, and \code{known.helix}.
	}
	\item{arcs}{
		TRUE if the structure should be plotted as arcs.  Arcs may be styled
		with styling columns, see example and \code{\link{plotHelix}} for
		details.
	}
	\item{add}{
		TRUE if graphical elements are to be added to an existing device, else
		a new plotting device is created with \code{\link{blankPlot}}.
	}
	\item{grid}{
		TRUE if the multiple sequence alignment is to be drawn as a grid of
		bases, else the multiple sequence alignment is drawn as equidistant
		horizontal lines.
	}
	\item{text}{
		Only applicable when grid is TRUE.  TRUE if the grid is to be filled
		with nucleotide character.
	}
	\item{legend}{
		TRUE if legend are to be shown.
	}
	\item{species}{
		If a number greater than 0 is given, then species names for the multiple
		sequence alignment will be printed along the left side.  This name is
		typically the entire header lines of FASTA entries from
		\code{\link{readFasta}}, and can be manually manipulated using the
		\code{\link{names}} function.  The number specifies the start position
		relative to the left edge of the multiple sequence alignment).
	}
	\item{base.colour}{
		TRUE if bases are to be coloured by nucleotide instead of basepair
		conservation.
	}
	\item{palette}{
		A list of colour names to override the default colour palette.  When
		base.colour is TRUE, the first 6 colours will be used for colouring
		bases A, U, G, C, - (gap), and ? (everything else), respectively.  When
		base.colour is FALSE, the first 7 colours will be used for colouring
		conserved basepairs, covarying basepairs, one-sided conserved basepairs,
		invalid basepairs, unpaired bases, gaps, and bases/pairs with ambiguous
		bases, resepctively.  If the palette is shorter than the expected
		length, the palette will simply cycle.  \dQuote{NA} is a valid colour,
		that will effectively plot nothing.
	}
	\item{flip}{
		If TRUE, the entire plot will be flipped upside down.  Note that this
		is \emph{not} a perfect mirror image about the horizon.
	}
	\item{grid.col, grid.lwd}{
		The colour and line width of the borders displayed when \code{grid} is
		TRUE.
	}	
	\item{text.cex, text.col, text.font, text.family}{
		cex, col, family and font for the text displayed via the \code{text}
		option.  Use \code{help("par")} for more information the paramters.
	}
	\item{species.cex, species.col, species.font, species.family}{
		cex, col, family and font for the species text displayed via the
		\code{species} option.  Use \code{help("par")} for more information the
		paramters.
	}
	\item{shape}{
		One of "circle", "triangle", or "square", specifying the shape of
		the arcs.
	}
	\item{y}{
		The vertical offset of the entire figure relative to 0.
	}
	\item{conflict.lty, conflict.col}{
		Determines the line type (style) and colour to be used for conflicting
		basepairs.  By default, conflicting helices are drawn as dotted lines
		(\code{lty = 2}) and whatever colour was originally assigned to it
		(\code{col = NA}). Conflicting helices may be coloured by setting
		\code{conflict.col} to some R-compatible colour name.  If both arguments
		are set to \code{NA}, then no attempt to exclude conflicting helices 
		will be made when colouring covariance plot columns, which in most cases
		will render the plot nonsensical.
	}
	\item{miss}{
		The colour for unpredicted arcs in overlapping diagrams, see
		\code{plotOverlapHelix} for more information.
	}
	\item{any}{
		Passed to \code{\link{plotOverlapHelix}}.  Determines definition of
		\dQuote{overlap}, if TRUE, then any basepair overlapping makes a helix
		overlap, if FALSE, then all basepairs in a helix must overlap to count
		as overlapping.
	}
	\item{pad}{
		A four integer array passed to \code{\link{blankPlot}}, specifies the
		number of pixels to pad the bottom, left, top and right sides of the
		figure with, repsectively.
	}
	\item{\dots}{
		In \code{plotCovariance}, these are additional arguments passed to
		\code{\link{blankPlot}}, useful arguments include \sQuote{lwd},
		\sQuote{col}, \sQuote{cex} for line width, line colour, and text size,
		respectively. \code{help('par')} for more.

		For \code{plotDoubleCovariance} and \code{plotOverlapCovariance}, these
		are additional arguments passed to \code{plotCovariance} (and thus
		indirectly also to \code{\link{blankPlot}}).
	}
}
\seealso{
	\code{\link{plotHelix}}

	\code{\link{plotDoubleHelix}}

	\code{\link{plotOverlapHelix}}

	\code{\link{colourByCovariation}}

	\code{\link{colourByConservation}}

	\code{\link{colourByCanonical}}
}
\author{
	Daniel Lai
}
\examples{
	data(helix)

	# Basic covariance plot
	plotCovariance(fasta, known, cex = 0.8, lwd = 1.5)

	# Grid mode
	plotCovariance(fasta, known, grid = TRUE, text = FALSE, cex = 0.8)

	# Global style and nucleotide colouring
	plotCovariance(fasta, known, grid = TRUE, text = FALSE, base.colour = TRUE)

	# Styling indivual helices with styling columns
	known$col <- c("red", "blue")
	plotCovariance(fasta, known, lwd = 2, cex = 0.8)

	# Use in combination with colourBy functions
	cov <- colourByCovariation(known, fasta, get = TRUE)
	plotCovariance(fasta, cov)
	legend("topleft", legend = attr(cov, "legend"),
		fill = attr(cov, "fill"), title = "Covariation")

}
\keyword{aplot}

