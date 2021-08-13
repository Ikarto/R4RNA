\name{internal}
\title{R4RNA Package Internal Functions}
\alias{R4RNA-internal}
\alias{openFileOrConnection}
\alias{anglesArc}
\alias{arc}
\alias{getBaseColours}
\alias{getCovarianceColours}
\alias{getSequenceColour}
\alias{plotCovarianceLine}
\alias{plotArrow}
\alias{expandNumberedHelix}
\alias{readMSA}

\description{
	Use only if you know what you're doing :)
	
	Accessible via R4RNA:::function_name
}

\usage{
	% IO helper functions
	openFileOrConnection(file, mode = "r")

	% Plotting helper functions
	anglesArc(v, theta) 
	arc(c, r, v, theta, ...)
	getBaseColours(msa)
	getCovarianceColours(msa, helix)
	getSequenceColour(msa, helix)
	plotCovarianceLine(x, y, cols, ...)
	plotArrow(start, y = 0, length = start * 0.03)

	% Manip helper functions
	expandNumberedHelix(helix)
	
}
\author{
	Daniel Lai
}

\keyword{internal}

