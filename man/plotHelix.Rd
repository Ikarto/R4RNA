\name{Plot Helix Structures}
\alias{plotHelix}
\alias{plotDoubleHelix}
\alias{plotOverlapHelix}
\alias{plotArcs}
\alias{plotArc}
\title{
	Plots helices in arc diagram
}
\description{
	Plots a helix data frame as an arc diagram, with styling possible with
	properly named additional columns on the data frame. 
}
\usage{
	plotHelix(helix, y = 0, flip = FALSE, line = FALSE, arrow = FALSE,
		add = FALSE, shape = "circle", ...)
	plotDoubleHelix(top, bot, line = TRUE, arrow = FALSE, add = FALSE, ...)
	plotOverlapHelix(predict, known, miss = "black", line = TRUE,
		arrow = FALSE, add = FALSE, any = FALSE, ...)
	plotArcs(i, j, length, y = 0, flip = FALSE, shape = "circle", ...)
	plotArc(i, j, y = 0, flip = FALSE, shape = "circle", ...)

}
\arguments{
	\item{helix, top, bot, predict, known}{
		Helix data.frames, with the four mandatory columns.  Any other column
		will be considered a styling column, and will be used for styling the
		helix.  See \code{example} for styling usage.  See \code{Details} for
		exact usage of each helix.
	}
	\item{y}{
		The vertical offset of the arc base relative to 0 along the y-axis.
	}
	\item{flip}{
		If TRUE, flips the arcs upside down about the y-axis.
	}
	\item{line}{
		If TRUE, a horizontal line representing the sequence is plotted.
	}
	\item{arrow}{
		If TRUE, an arrow is played on the right end of the line.
	}
	\item{add}{
		If TRUE, graphical elements are added to the active plot device, else
		a new plot device is created for the plot.
	}
	\item{shape}{
		One of "circle", "triangle", or "square", specifying the shape of
		the arcs.
	}
	\item{miss}{
		The colour for unpredicted arcs in overlapping diagrams, see
		\code{details} for more information.
	}
	\item{any}{
		Passed to \code{\link{isOverlappingHelix}}.  Determines definition of
		\dQuote{overlap}, if TRUE, then any basepair overlapping makes a helix
		overlap, if FALSE, then all basepairs in a helix must overlap to count
		as overlapping.
	}
	\item{i, j}{
		The starting and ending position of the arc along the x-axis
	}
	\item{length}{
		The total number of arcs to draw by incrementing \emph{i} and
		decrementing \emph{j}.  Used to draw helices.
	}
	\item{\dots}{
		Any additional parameters passed to \code{par}
	}
}
\details{
  \code{plotHelix} creates a arc diagram with all arcs on top,
  \code{plotDoubleHelix} creates a diagram with arcs on the top and bottom.
  \code{plotOverlapHelix} is slight trickier, and given two structures
  \code{predict} and \code{known}, plots the predicted helices that are known
  on top, predicted helices that are not known on the bottom, and finally plots
  unpredicted helices on top in the colour defined by \code{miss}.
  
  \code{plotArc} and \code{plotArcs} are the core functions that make everything
  work, and may be used for extreme fine-tuning and customization.
}

\seealso{
	\code{\link{colourByCount}}
}
\author{
	Daniel Lai
}
\examples{
	data(helix)

	# Plot helix plain
	plotHelix(known)

	# Apply global appearance options
	plotHelix(known, line = TRUE, arrow = TRUE, col = "blue", lwd = 1.5)

	# Add extra column with styling options
	known$lty <- 1:4
	known$lwd <- 1:2
	known$col <- c(rgb(1, 0, 0), "orange", "yellow", "#00FF00", 4, "purple")
	plotHelix(known)

	# Manually colour helices according to value
	helix$col <- "red"
	helix$col[which(helix$value < 1e-3)] <- "orange"
	helix$col[which(helix$value < 1e-4)] <- "green"
	helix$col[which(helix$value < 1e-5)] <- "blue"
	plotHelix(helix)

	# Automatically creating a similar plot with legend
	coloured <- colourByValue(helix, log = TRUE, get = TRUE)
	plotHelix(coloured, line = TRUE, arrow = TRUE)
	legend("topleft", legend = attr(coloured, "legend"),
		fill = attr(coloured, "fill"), title = "P-value", text.col = "black")	

	# Plot both helices with styles
	plotDoubleHelix(helix, known)

	# Overlap helix
	plotOverlapHelix(helix, known)
}
\keyword{aplot}
