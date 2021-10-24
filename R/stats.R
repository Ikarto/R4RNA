## Copyright (C) 2011 Daniel Lai and Irmtraud M. Meyer (www.e-rna.org)
## Contact information: irmtraud.meyer@cantab.net

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @export
#'
colourByCovariation <- function(helix, msa, cols, get = FALSE) {
	output <- colourBy.internal(helixCovariation(helix, msa), 2, -2, cols)
	if (get) {
		helix$col <- output
		attr(helix, "legend") <- attr(output, "legend")
		attr(helix, "fill") <- attr(output, "fill")
		return(helix)
	} else {
		return(output)
	}
}


#' @export
#'
colourByConservation <- function(helix, msa, cols, get = FALSE) {
	output <- colourBy.internal(helixConservation(helix, msa), 1, 0, cols)
	if (get) {
		helix$col <- output
		attr(helix, "legend") <- attr(output, "legend")
		attr(helix, "fill") <- attr(output, "fill")
		return(helix)
	} else {
		return(output)
	}
}



#' @export
#'
colourByCanonical <- function(helix, msa, cols, get = FALSE) {
	output <- colourBy.internal(helixCanonical(helix, msa), 1, 0, cols)
	if (get) {
		helix$col <- output
		attr(helix, "legend") <- attr(output, "legend")
		attr(helix, "fill") <- attr(output, "fill")
		return(helix)
	} else {
		return(output)
	}
}



#' @export
#'
colourBy.internal <- function(values, start, end, cols) {
	if (missing(cols)) {
		cols <- defaultPalette()
	}
	cols <- rev(cols)
	seq <- seq(start, end, length.out = length(cols) + 1)
	levels <- cut(values, breaks = seq, include.lowest = TRUE)
	legend <- levels(levels)
	fill <- cols
	output <- cols[as.numeric(levels)]
	attr(output, "legend") <- legend
	attr(output, "fill") <- fill
	return(output)
}



#' @export
#'
# TODO: Three following functions could probably be refactored?
# Given a helix data frame x, returns a list of covariation values corresponding to each row
helixCovariation <- function(helix, msa) {
	if (!is.helix(helix)) {
		stop("Invalid input")
	}
	ex <- expandHelix(helix)
	num <- rep(1:nrow(helix), times = helix$length)
	cov <- c()
	for (i in 1:nrow(ex)) {
		cov[i] <- basepairCovariation(msa, ex$i[i], ex$j[i])
	}
	pt <- 1
	out = c()
	for (i in 1:nrow(helix)) {
		out[i] <- mean(cov[seq(pt, length.out = helix$length[i])])
		pt <- pt + helix$length[i]
	}
	return(out)
}



#' @export
#'
# Given a helix data frame x, returns a list of covariation values corresponding to each row
helixConservation <- function(helix, msa) {
	if (!is.helix(helix)) {
		stop("Invalid input")
	}
	ex <- expandHelix(helix)
	num <- rep(1:nrow(helix), times = helix$length)
	cons <- c()
	for (i in 1:nrow(ex)) {
		cons[i] <- basepairConservation(msa, ex$i[i], ex$j[i])
	}
	pt <- 1
	out = c()
	for (i in 1:nrow(helix)) {
		out[i] <- mean(cons[seq(pt, length.out = helix$length[i])])
		pt <- pt + helix$length[i]
	}
	return(out)
}



#' @export
#'
helixCanonical <- function(helix, msa) {
	if (!is.helix(helix)) {
		stop("Invalid input")
	}
	ex <- expandHelix(helix)
	num <- rep(1:nrow(helix), times = helix$length)
	cano <- c()
	for (i in 1:nrow(ex)) {
		cano[i] <- basepairCanonical(msa, ex$i[i], ex$j[i])
	}
	pt <- 1

	out = c()
	for (i in 1:nrow(helix)) {
		out[i] <- mean(cano[seq(pt, length.out = helix$length[i])])
		pt <- pt + helix$length[i]
	}
	return(out)
}



#' @export
#'
# Given a multiple sequence alignment, and two column indices indicating the 5' partner, and 3' partner respectively,
# calculates the covariation of a base pair at this position.
# values ranges between -2 and 2
basepairCovariation <- function(msa, pos.5p, pos.3p) {
	covariation <- 0
	seq.pair.count <- choose(length(msa), 2)
	base.pairs <- paste(substr(msa, pos.5p, pos.5p), substr(msa, pos.3p, pos.3p), sep="")

	pair.counts <- table(base.pairs)
	pair.names <- names(pair.counts)

	if (length(pair.counts) <= 1) {
		return(0)
	}

	for (pairA in 2:length(pair.counts)) {
		for (pairB in 1:(pairA-1)) {
			bpA <- pair.names[pairA]
			bpB <- pair.names[pairB]
			occurences <- pair.counts[pairA] * pair.counts[pairB]
			hdist <- hamming(bpA, bpB)
			coef <- ifelse(isValidBp(bpA) & isValidBp(bpB), 1, -1)
			covariation <- covariation + (coef * hdist * occurences)
		}
	}
	return (covariation / seq.pair.count)
}

# returns the hamming distance between strings s1 and s2
# hamming distance = number of positions which differ between s1 and s2
# ASSUMES same length
hamming <- function(s1, s2) {
	return(sum(strsplit(s1, "")[[1]] != strsplit(s2, "")[[1]]))
}

isValidBp <- function(bp) {
	return(bp %in% c("GC", "CG", "UA", "AU", "GU", "UG", "TA", "AT", "TG", "GT"))
}



#' @export
#'
basepairConservation <- function(msa, pos.5p, pos.3p) {
	return((baseConservation(msa, pos.5p) + baseConservation(msa, pos.3p)) / 2)
}


#' @export
#'
baseConservation <- function(msa, pos) {
	return(columnPercentIdentity(substr(msa, pos, pos)))
}



#' @export
#'
basepairCanonical <- function(msa, pos.5p, pos.3p) {
	base.pairs <- paste(substr(msa, pos.5p, pos.5p), substr(msa, pos.3p, pos.3p), sep="")
	canonical.count <- 0
	for (bp in base.pairs) {
		if (isValidBp(bp)) {
			canonical.count <- canonical.count + 1
		}
	}
	return(canonical.count/length(base.pairs))
}




#' @export
#'
# linear implementation of column percent identity
# does not work... cannot correctly apply denominator
columnPercentIdentity <- function(bases) {
	identities <- 0
	seq.pair.count <- choose(length(bases), 2)

	base.counts <- table(bases)
	base.names <- names(base.counts)

	if (sum(base.counts) != length(bases)) {
		warning("columnPercentIdentity: length of column does not equal sum(base count)")
	}

	num.gap.pairs <- 0
	for (base.idx in 1:length(base.counts)) {
		count <- base.counts[base.idx]
		base <- base.names[base.idx]
		if (base == "-") {
			num.gap.pairs <- choose(count, 2)
		} else {
			identities <- identities + choose(count, 2)
		}
	}
	return(identities / (seq.pair.count - num.gap.pairs))
}



#' @export
#'
# given a multiple sequence alignment, and a helix data structure,
# this function returns the covariation of the entire alignment.
# the formula used is described on pg. 86 of Nick Wiebe's MSc thesis
alignmentCovariation <- function(msa, helix) {
	if (any(helix$length != 1)) {
		helix <- expandHelix(helix)
	}
	covariation <- sum(apply(helix, 1, function(bp, msa) {
			basepairCovariation(msa, bp['i'], bp['j'])
		}, msa))
	return (covariation/nrow(helix))
}



#' @export
#'
# calculates average percent identity of all possible pairs of sequences in the alignment
# percent identity is equal to number of non-gap matches / number of non-gap columns
# definition of percent identity: average over all seq pairs (# identities / (aligned positions + internal gaps))
alignmentConservation <- function(msa) {
	seq.pair.count <- choose(length(msa), 2)
	pid.sum <- 0
	msa.bases <- strsplit(msa, "")
	for (seqB.index in 2:length(msa)) {
		for (seqA.index in 1:(seqB.index-1)) {
			seqA.bases <- msa.bases[[seqA.index]]
			seqB.bases <- msa.bases[[seqB.index]]
			alignedPositions <- which(seqA.bases != "-" & seqB.bases != "-")
			start <- min(alignedPositions)
			end <- max(alignedPositions)
			seqA.bases <- seqA.bases[start:end]
			seqB.bases <- seqB.bases[start:end]
			contains.nucleotide <- which(seqA.bases != "-" | seqB.bases != "-")
			seqA.cleaned <- seqA.bases[contains.nucleotide]
			seqB.cleaned <- seqB.bases[contains.nucleotide]
			identities <- sum(seqA.cleaned == seqB.cleaned)
			pid.sum <- pid.sum + identities/length(seqA.cleaned)
		}
	}
	return (pid.sum/seq.pair.count)
}



#' @export
#'
# calculates percent of base pairs that are canonical (GC, CG, AU, UA, GU, UG) in the alignment
alignmentCanonical <- function(msa, helix) {
	if (any(helix$length != 1)) {
		helix <- expandHelix(helix)
	}
	percentCanonical <- 0
	for (bp.num in 1:nrow(helix)) {
		bp.i <- helix[bp.num, 'i']
		bp.j <- helix[bp.num, 'j']
		percentCanonical <- percentCanonical +
			basepairCanonical(msa, bp.i, bp.j)
	}
	return(percentCanonical / nrow(helix))
}


#' @export
#'
structureMismatchScore <- function(msa, helix) {
	if (!is.helix(helix)) {
		stop("Invalid helix data.frame")
	} else {
		helix <- expandHelix(helix)
	}
	return(unlist(lapply(msa, structureMismatchScoreInternal, helix)))
}


structureMismatchScoreInternal <- function(sequence, helix) {
	base.5p <- substring(sequence, helix$i, helix$i)
	base.3p <- substring(sequence, helix$j, helix$j)
	gaps.5p <- base.5p == "-"
	gaps.3p <- base.3p == "-"
	one.sided.gap <- xor(gaps.5p, gaps.3p)
	two.sided.gap <- gaps.5p & gaps.3p
	invalid.pair <- !isValidBp(paste(base.5p, base.3p, sep = ""))
	non.canonical.pair <- invalid.pair & !(one.sided.gap | two.sided.gap)
	return(c(2 * length(which(one.sided.gap)) + length(which(two.sided.gap))
		+ length(which(non.canonical.pair))))
}



#' @export
#'
alignmentPercentGaps <- function(msa) {
	return(unlist(lapply(lapply(strsplit(msa, ""), '==', '-'), sum)) / nchar(msa))
}

