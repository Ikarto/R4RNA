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
as.helix <- function(x, length = NULL) {
	x <- as.data.frame(x)
	if (nrow(x) == 0) {
		x <- data.frame(i = vector(), j = vector(), length = vector(),
			value = vector())
	}

	if(ncol(x) < 4) {
		stop("Expecting data.frame-like structure with at least 4 columns")
	}
	names <- colnames(x)
	names[1:4] <- c("i", "j", "length", "value")
	colnames(x) <- names
	rownames(x) <- NULL
	x$i <- as.integer(x$i)
	x$j <- as.integer(x$j)
	x$length <- as.integer(x$length)
	x$value <- as.numeric(x$value)

	if (is.null(length)) {
		length <- max(x$i, x$j)
	}

	attr(x, "length") <- length
	return(x)
}


#' @export
#'
is.helix <- function(x) {
	valid <- TRUE
	if (!is.data.frame(x)) {
		warning("Not a valid data frame")
		valid <- FALSE
	}
	if (ncol(x) < 4) {
		warning("Should have at least 4 columns")
		valid <- FALSE
	}
	if (is.null(attr(x, "length"))) {
		warning("No length attribute")
		valid <- FALSE
	}
	if (!all(lapply(x, class)[1:4] == c("integer", "integer", "integer",
			"numeric"))) {
		warning("Columns of invalid class")
		valid <- FALSE
	}
	if(!all(colnames(x)[1:4] == c("i", "j", "length", "value"))) {
		warning("Invalid column names")
		valid <- FALSE
	}
	return(valid)
}


#' @export
#'
collapseHelix <- function(helix) {
	if(!is.helix(helix)) {
		stop("Input is not a valid helix, aborting")
	} else {
		if (nrow(helix) == 0) {
			return(helix)
		}
		helix <- expandHelix(helix)
		helix <- helix[!duplicated(helix), ]
	}
	output <- data.frame()
	length <- attr(helix, "length")
	values <- unique(helix$value)
	for(v in 1:length(values)) {
		value <- values[v]
		if (is.na(value)) {
			y <- helix[which(is.na(helix$value)), ]
		} else {
			y <- helix[which(helix$value == value), ]
		}
		y$sum <- y$i + y$j
		y <- y[order(y$sum, y$i), ]

		count <- nrow(y)
		if (y$i[1] == y$i[count] - (count - 1) &&
			y$j[1] == y$j[count] + (count - 1)) {
			output <- rbind(output, c(y$i[1], y$j[1], count, value))
		} else {
			y$start <- 0
			y$start[1] <- 1
			expecti <- y$i[1]
			expectj <- y$j[1]
			stem <- c(expecti, expectj, 0)
			for (i in 1:count) {
				if (y$i[i] == expecti && y$j[i] == expectj) {
					expecti <- expecti + 1
					expectj <- expectj - 1
					stem[3] <- stem[3] + 1
				} else {
					output <- rbind(output, c(stem, value))
					y$start[i] <- 1
					stem <- c(y$i[i], y$j[i], 1)
					expecti <- y$i[i] + 1
					expectj <- y$j[i] - 1
				}
			}
			output <- rbind(output, c(stem, value))
		}
	}
	return(as.helix(output, length))
}


#' @export
#'
expandHelix <- function(helix) {
	if(!is.helix(helix)) {
		stop("Input is not a valid helix, aborting")
	} else {
		if (nrow(helix) == 0 || all(helix$length == 1)) {
			return(helix)
		}
	}
	length <- attr(helix, "length")
	output <- data.frame(i = 0, j = 0, length = 0, value = 0)
	itr <- 1
	for(i in 1:nrow(helix)) {
		start <- helix$i[i]
		end <- helix$j[i]
		width <- helix$length[i]
		value <- helix$value[i]
		for(offset in 0:(width - 1)) {
			output[itr, ] <- c(start + offset, end - offset, 1, value)
			itr <- itr + 1
		}
	}
	return(as.helix(output, length))
}


expandNumberedHelix <- function(helix) {
	number <- rep(1:nrow(helix), helix$length)
	helix <- expandHelix(helix)
	helix$number <- number
	return(helix)
}


#' @export
#'
isConflictingHelix <- function(helix, any = TRUE) {
	if (!is.helix(helix)) {
		stop("Not a valid helix data frame, aborting...")
	}
	helices <- nrow(helix)
	if(helices == 0) {
		return (vector());
	}

	helix <- expandNumberedHelix(helix)
	end <- nrow(helix)
	helix$keep <- T
	for (i in 1:(end - 1)) {
		if (helix$keep[i]) {
			kill <- which(helix$i[(i + 1):end] == helix$i[i] |
				helix$j[(i + 1):end] == helix$j[i] |
				helix$j[(i + 1):end] == helix$i[i] |
				helix$i[(i + 1):end] == helix$j[i])
			helix$keep[kill + i] <- F
		}
	}
	keep <- c()
	for (i in 1:helices) {
		if (any) {
			keep <- c(keep, all(helix$keep[which(helix$number == i)]))
		} else {
			keep <- c(keep, any(helix$keep[which(helix$number == i)]))
		}
	}
	return(!keep)
}


#' @export
#'
isDuplicatingHelix <- function(helix, any = TRUE) {
	if (!is.helix(helix)) {
		stop("Invalid helix data frame, aborting")
	}
	helices <- nrow(helix)

	if(helices == 0) {
		return (vector());
	}

	helix <- expandNumberedHelix(helix)
	flip <- which(helix$i > helix$j)
	tmp <- helix$i[flip]
	helix$i[flip] <- helix$j[flip]
	helix$j[flip] <- tmp

	helix$keep <- !duplicated(helix[, c("i", "j")])
	keep <- c()
	for (i in 1:helices) {
		if (any) {
			keep <- c(keep, all(helix$keep[which(helix$number == i)]))
		} else {
			keep <- c(keep, any(helix$keep[which(helix$number == i)]))
		}
	}
	return(!keep)
}



#' @export
#'
isOverlappingHelix <- function(helix, query, any = FALSE) {
	if(!is.helix(helix) | !is.helix(query)) {
		stop("One of the inputs was not a valid helix data frame, aborting")
	}
	helices <- nrow(helix)

	if (helices == 0) {
		return (vector());
	}

	helix <- expandNumberedHelix(as.helix(helix[, 1:4]))
	y <- expandHelix(query)[, 1:4]
	y$length <- NULL
	y$value <- NULL
	y$hit <- T
	helix <- merge(helix, y, all.x = TRUE)
	helix$hit[is.na(helix$hit)] <- F
	overlap <- c()
	for (i in 1:helices) {
		if (any) {
			overlap <- c(overlap, any(helix$hit[which(helix$number == i)]))
		} else {
			overlap <- c(overlap, all(helix$hit[which(helix$number == i)]))
		}
	}
	return(overlap)
}



#' @export
#'
logseq <- function(from, to, length.out) {
	if (missing(length.out)) {
		length.out <- log10(max(from, to) / min(from, to)) + 1
	}
	return(exp(seq(log(from), log(to), length.out = length.out)))
}

logfloor <- function(x) {
	return(10 ** floor(log10(x)))
}

logceiling <- function(x) {
	return(10 ** ceiling(log10(x)))
}


#' @export
#'
basepairFrequency <- function(helix) {
	if (!is.helix(helix)) {
		stop("Not a valid helix data frame, aborting")
	}
	if (nrow(helix) == 0) {
		return(helix)
	} else {
		basepairs <- expandHelix(helix)[, c("i", "j")]
	}
	counts <- table(paste(basepairs$i, basepairs$j))
	pos <- as.integer(unlist(strsplit(names(counts), " ")))
	odds <- seq(1, length(pos), by = 2)
	output <- data.frame(i = pos[odds], j = pos[odds + 1], length = 1,
		value = as.integer(counts))
	output <- output[order(-output$value, output$i, output$j), ]
	return(as.helix(output, attr(helix, "length")))
}


# Returns a logical row each row of the matrix, TRUE of pseudoknotted, else false
# Assumes these are  helices of LENGTH 1
is_pseudoknotted <- function(row, matrix) {
	return((row$i <= matrix$i & matrix$i <= row$j & row$j <= matrix$j) |
	(matrix$i <= row$i & row$i <= matrix$j & matrix$j <= row$j))
}


#' @export
#'
# Adds a group column, in which helices of each numbered group are non-pseudoknotted
unknottedGroups <- function(helix) {
	if (!is.helix(helix)) {
		stop("Invalid input")
	}
	if (any(helix$length > 1)) {
		warning("Expanding helix to basepairs")
		helix <- expandHelix(helix)
	}
	if (nrow(helix) == 0) {
		stop("No basepairs in helix")
	}
	group <- rep(0, nrow(helix));
 	group[1] <- 1
	if (nrow(helix) == 1) {
		return(group)
	}
	for (i in 2:nrow(helix)) {
		test <- 1
		while (any(is_pseudoknotted(helix[i, ], helix[which(group == test), ]))) {
			test <- test + 1
		}
		group[i] <- test
	}
	return(group)
}

