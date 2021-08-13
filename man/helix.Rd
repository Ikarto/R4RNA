\name{Example Data}
\docType{data}
\alias{helix}
\alias{fasta}
\alias{known}
\title{Helices predicted by TRANSAT with p-values}
\description{
	This data set contains two sets of helices and a multiple sequence alignment. 
	The two sets of helices are \code{helices} and \code{known} which are
	helices predicted to occur for RNA sequence RF00458 by the program TRANSAT,
	and experimentally proposed structure of the same sequence, respectively.
	\code{fasta} is the seed homologues for the multiple sequence alignment
	obtained from the RFAM database.
}
\usage{
	data(helix)
}
\format{
	\code{helix} and \code{known} are 4 column data frames, where columns
	\bold{i} and \bold{j} denote the left-most and right-most basepairs, the
	\bold{length} is the number of \emph{consecutive} basepairs the helix
	contains, and the \bold{value} is assigned to each helix on a row.
	
	\code{fasta} is an array of named characters of length 7.
}
\references{
	Wiebe NJ, Meyer IM. (2010) \emph{TRANSAT-- method for detecting the
	conserved helices of functional RNA structures, including transient,
	pseudo-knotted and alternative structures}.  PLoS Comput Biol.
	6(6):e1000823.

	Gardner PP, Daub J, Tate J, Moore BL, Osuch IH, Griffiths-Jones S, Finn RD,
	Nawrocki EP, Kolbe DL, Eddy SR, Bateman A. (2011) \emph{Rfam: Wikipedia,
	clans and the ``decimal'' release}.  Nucleic Acids Res. 39(Database
	issue):D141-5.
}
\keyword{datasets}
