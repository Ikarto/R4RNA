\name{Write Helix}
\alias{writeHelix}
\title{
	Write out a helix data frame into a text file
}
\description{
	Write out a helix data frame into a text file into the four-column
	tab-delimited format with proper header and column names.
}
\usage{
	writeHelix(helix, file = stdout())
}
\arguments{
	\item{helix}{
		A helix data frame.
	}
	\item{file}{
		A character string pointing to a file path, or a file connection.
		Defaults to the console.
	}
}
\author{
	Daniel Lai
}
\examples{
	# Create helix data frame
	helix <- data.frame(2, 8, 3, 0.5)
	helix[2, ] <- c(5, 15, 4, -0.5)
	helix <- as.helix(helix)

	writeHelix(helix)
}
\keyword{IO}
\keyword{file}
