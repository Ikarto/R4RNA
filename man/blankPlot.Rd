\name{Create Blank Plot}
\alias{blankPlot}
\alias{maxHeight}
\title{
	Create a blank plotting canvas
}
\description{
	Creates a blank plotting canvas with the given dimensions, along with
	functions to find best values for the canvas dimensions.
}
\usage{
	blankPlot(width, top, bottom, pad = c(0, 0, 0, 0), scale = TRUE,
		scale.lwd = 1, scale.col = "#DDDDDD", scale.cex = 1, debug = FALSE,
		png = NA, pdf = NA, factor = ifelse(!is.na(png), 8, 1/9),
		no.par = FALSE, ...)
	maxHeight(helix)
}
\arguments{
	\item{width}{
		A number indicating the horizontal width of the blank plot.
	}
	\item{top, bottom}{
		The maximum and minimum values vertically to be displayed in the plot.
	}
	\item{pad}{
		An array of 4 integers, specifying the pixels of whitespace to pad
		beyond the dimensions given by top, bottom, and width.  Four number
		corresponding to padding on the bottom, left, top and right,
		respectively.  Default is c(0, 0, 0, 0).
	}
	\item{scale}{
		If TRUE, inserts a scale on the plot.
	}
	\item{scale.lwd, scale.col, scale.cex}{
		Allows manual modification of the scale's line width and colour,
		respectively.
	}
	\item{png, pdf}{
		If one or the other is set to a filename, a file in png or pdf format
		will be produced respectively.  If both are set to non-NA values,
		\code{png} will have priority.
	}
	\item{factor}{
		The scaling factor used to produce plots of png or pdf format.  Should
		be set so after multiplication of the \code{top, bot}, etc arguments,
		good document dimension in pixels with \code{png} and inches for
		\code{pdf} will be produced.
	}
	\item{debug}{
		If TRUE, frames the boundaries of the intended plotting space in red,
		used to determine if inputs produce expected output area.  Also outputs
		to STDIN dimensions of the plot.
	}
	\item{no.par}{
		Suppresses the internal call to \code{par} in the function if set to
		TRUE, useful for using \code{par} arguments such as \code{mfrow}, etc.
	}
	\item{\dots}{
		Additional arguments passed to \code{\link{par}} when \code{no.par}
		is FALSE, common ones include \sQuote{lwd}, \sQuote{col}, \sQuote{cex}
		for line width, line colour, and text size, respectively.
		\code{help('par')} for more.  When \code{no.par} is set to TRUE, this
		option does nothing, and manually calling \code{par} is required prior
		to the calling of this function.
	}
	\item{helix}{
		A helix data.frame
	}
}
\details{
	\code{blankPlot} creates a blank plot with the given dimensions, with
	minimal margins around the plot and no axis or labels.  If more control
	is required, using \code{\link{plot}} directly would be more efficient.

	\code{maxHeight} returns the height that the highest helix would require,
	and can be used to determine \code{top} and \code{bottom} for
	\code{blankPlot}.
}
\value{
	\code{maxHeight} returns a numeric integer.
}
\seealso{
	\code{\link{plotHelix}}
}
\author{
	Daniel Lai
}

\examples{
	# Create helix and obtain height
	helix <- as.helix(data.frame(1, 37, 12, 0.5))
	height <- maxHeight(helix)
	print(height)

	# Use height to create properly sized plot
	width <- attr(helix, "length")
	blankPlot(width, height, 0)

	# Add helix to plot
	plotHelix(helix, add = TRUE)
}
\keyword{aplot}

