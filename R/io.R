## Copyright (C) 2020-2021 Volodymyr Tsybulskyi and Irmtraud M. Meyer (www.e-rna.org)
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
#'
writeHelix <- function(helix, file = stdout()) {
	if (!is.helix(helix)) {
		warning("Data frame not a compliant helix, reformatting...")
		helix <- as.helix(helix, attr(helix, "length"))
	}
	writeLines(paste("#", attr(helix, "length")), file)
	suppressWarnings(write.table(helix, file, quote = F, row.names = F,
		sep = "\t", append = T))
}


#' @export
#'
writeFasta <- function(msa, file = stdout(), wrap = NA) {
	desc <- names(msa)
	if (!is.na(wrap)) {
		if (wrap < 1) {
			stop("Wrap must be a non-zero integer")
		}
		if (length(unique(nchar(msa))) == 1) {
			length <- base:nchar(msa[1])
			start <- seq(1, length, by = wrap)
			end <- start + wrap - 1
			split <- lapply(msa, substring, start, end)
			msa <- unlist(lapply(split, paste, collapse = "\n"))
		} else {
			start <- lapply(nchar(msa), seq, from = 1, by = wrap)
			end <- lapply(start, '+', wrap - 1)
			for (i in 1:length(msa)) {
				split <- substring(msa[i], start[[i]], end[[i]])
				msa[i] <- paste(split, collapse = "\n")
			}
		}
	}
	text <- paste(">", desc, "\n", msa, sep = "")
	writeLines(text, file)
}


#' @export
#'
readConnect <- function(file) {
	con <- openFileOrConnection(file)
	on.exit(close(con))
	output <- data.frame()
	max <- 0
	while(length(header <- readLines(con, n = 1)) > 0) {
		width <- strsplit(header, "\\s+")[[1]]
		width <- width[which(width != "")]
		if (length(width) == 0) { next }
		width <- suppressWarnings(as.numeric(width[1]))
		if (is.na(width)) {
			stop("Structure length missing in header of connect file")
		}
		energy <- regexpr("(ENERGY|dG)\\s*=?\\s*\\S*", header, ignore.case = T)
		energy <- substr(header, energy, energy + attr(energy, "match.length") - 1)
		energy <- sub("(ENERGY|dG)\\s*=?\\s*", "", energy, ignore.case = T)
		energy <- suppressWarnings(as.numeric(energy))
		table <- read.delim(con, header = FALSE, nrows = width, sep = "")
		if (ncol(table) != 6) {
			print(ncol(table))
			print(table)
			stop("Requires 6 data columns in connect format, aborting.")
		}
		table <- table[which(table[, 5] != 0), c(1, 5)]
		names(table) <- c("i", "j")
		table$length <- 1L
		table$value <- energy
		table <- table[which(table$j > table$i), ]
		max <- max(width, max, na.rm = TRUE)
		output <- rbind(output, table)
	}
	return(as.helix(output, max))
}



#' @export
#'
readBpseq <- function(file) {
	con <- openFileOrConnection(file)
	on.exit(close(con))
	lines <- readLines(con)

	file_line <- lines[grep("^Filename", lines, ignore.case = TRUE)]
	organism_line <- lines[grep("^Organism", lines, ignore.case = TRUE)]
	accession_line <- lines[grep("^Accession", lines, ignore.case = TRUE)]
	bps <- lines[grep("^\\d+", lines)]
	cells <- strsplit(bps, "\\s+")

	if (unique(unlist(lapply(cells, length))) != 3) {
		stop("Expecting three column data in BPSEQ format")
	}

	matrix <- data.frame(matrix(unlist(cells), ncol = 3, byrow = T))
	i <- as.integer(as.character(matrix[, 1]))
	j <- as.integer(as.character(matrix[, 3]))
	width <- max(i, na.rm = TRUE)
	seq <- paste(matrix[1:width, 2], collapse = "")
	filename <- unlist(strsplit(file_line[1], ":\\s*"))[2]
	organism <- unlist(strsplit(organism_line[1], ":\\s*"))[2]
	accession <- unlist(strsplit(accession_line[1], ":\\s*"))[2]

	subset <- which(i < j)
	helix <- data.frame(i = i[subset], j = j[subset], length = 1L, value = NA)

	attr(helix, "organism") <- organism;
	attr(helix, "sequence") <- seq;
	attr(helix, "accession") <- accession;
	attr(helix, "filename") <- filename;
	return(as.helix(helix, width));
}



#' @export
#'
readHelix <- function(file) {
	file <- openFileOrConnection(file)
	on.exit(close(file))

	width <- strsplit(sub("#(\\s*)", "", readLines(file, n = 1)), "\\s")[[1]][1]
	width <- suppressWarnings(as.numeric(width))

	input <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
	colnames(input)[1:4] <- c("i", "j", "length", "value")

	max <- max(input$i, input$j, na.rm = TRUE)
	if (is.na(width)) {
		warning("Could not read valid length, using maximum position instead.")
		width <- max
	}
	if (width < max) {
		warning("Changing length to maximum basepair position.")
	}

	return(as.helix(input, width))
}


#' @export
#'
readVienna <- function(file, palette = NA) {
	if (length(grep("^>", readLines(file, n = 1))) == 1) {
		lines <- readFasta(file)
	} else {
		lines <- readLines(file)
	}
	name <- ""
	if (length(names(lines)) > 0) {
		name <- names(lines)[1];
	}

	raw <- lines[grep('[(].*[)]|<.*>|[[].*[]]|[{].*[}]', lines)]
	split <- strsplit(raw, "\\s")
	struct <- unlist(lapply(split, "[", 1))
	val <- unlist(lapply(split, "[", 2))
	seq <- ""

	length <- nchar(struct)
	halve <- grep("G|U", struct, ignore.case = TRUE)
	if (length(halve) > 0) {
		seq <- substr(struct[halve], 1, length/2);
	}

	struct[halve] <- substr(struct[halve], (length / 2) + 1, length)
	length <- nchar(struct)

	output <- data.frame()
	for (i in 1:length(struct)) {
		helix <- viennaToHelix(struct[i], val[i], palette = palette)
		output <- rbind(output, helix)
	}
	attr(output, "sequence") <- seq;
	attr(output, "name") <- name;
	return(as.helix(output, max(length)))
}



#' @export
#'
readFasta <- function (file, filter = FALSE) {
	file <- openFileOrConnection(file)
	on.exit(close(file))

    s1 <- scan(file = file, what = "", sep = "\n", quote = "",
        allowEscapes = FALSE, quiet = TRUE, comment.char = ";")
    descriptions <- which(substr(s1, 1L, 1L) == ">")
    numF <- length(descriptions)
    if (numF == 0)
        stop("no FASTA sequences found")
    dp <- descriptions + 1L
    dm <- descriptions - 1L
    end <- c(dm[-1], length(s1))
    tmp <- lapply(seq_len(numF), function(i) {
        desc <- s1[descriptions[i]]
        desc <- substr(desc, 2L, nchar(desc))
        if (end[i] >= dp[i]) {
            seq <- paste(s1[dp[i]:end[i]], collapse = "")
        }
        else {
            warning("record \"", desc, "\" contains no sequence")
            seq <- ""
        }
        list(desc = desc, seq = seq)
    })
    output <- as.character(unlist(lapply(tmp, "[[", "seq")))
    names(output) <- unlist(lapply(tmp, "[[", "desc"))

	if(filter) {
		output <- toupper(output)
		output <- gsub("T", "U", output)
		output <- gsub("[.]", "-", output)
	}
    return(output)
}



#' @export
#'
readMSA <- function(file, filter = TRUE) {
	message("readMSA it deprecated and will be removed, use readFasta(..., filter = TRUE)")
	return(readFasta(file, TRUE))
}


openFileOrConnection <- function(file, mode = "r") {
    if (is.character(file)) {
        file <- file(file, mode)
    }
    else {
        if (!inherits(file, "connection"))
            stop("'file' must be a character string or connection")
        if (!isOpen(file)) {
            open(file, mode)
        }
    }
    return(file)
}


#' @export
#'
viennaToHelix <- function(vienna, value = NA, palette = NA) {
	test <- regexpr("\\-?(\\d+\\.?\\d*|\\d*\\.?\\d+)", value)
	if (!is.na(test) & test != -1) {
		value <- as.numeric(substr(value, test,
			test + attr(test, "match.length") - 1))
	} else {
		value <- NA
	}
	brackets <- c('(' = ')', '<' = '>', '[' = ']', '{' = '}',
		'A' = 'a', 'B' = 'b', 'C' = 'c', 'D' = 'd')
	stacks <- list(')' = c(), '>' = c(), ']' = c(), '}' = c(),
		'a' = c(), 'b' = c(), 'c' = c(), 'd' = c())
	pairs <- c()
	valids <- c('(', ')', '<', '>', '[', ']', '{', '}', 'A',
		'a', 'B', 'b', 'C', 'c', 'D', 'd', '.', '-')
	chars <- unlist(strsplit(vienna, ""))
	chars <- chars[which(chars %in% valids)]
	chars[which(chars == '-')] <- '.';
	if (nchar(vienna) != length(chars)) {
		warning("Removed invalid character from input string")
	}

	output <- data.frame()
	for (j in 1:length(chars)) {
		char <- chars[j]
		if (char != '.') {
			if (!is.na(brackets[char])) {
				stacks[[brackets[char]]] <- append(stacks[[brackets[char]]], j)
			} else {
				if (!is.na(stacks[char]) && length(stacks[char]) > 0) {
					end <- length(stacks[[char]])
					partner <- stacks[[char]][end]
					if (length(partner) == 0) {
						stop(paste("Unbalanced", char,
							"bracket at position", j))
					}
					output <- rbind(output, c(partner, j, 1, value,
						which(brackets == char)))
					stacks[[char]] <- stacks[[char]][-end]
				} else {
					stop(paste("Unrecognized symbol", char,
						"in position", j))
				}
			}
		}
	}
	for (i in 1:length(brackets)) {
		if (length(stacks[[brackets[i]]]) > 0) {
			stop(paste("Unable to find closing partners for all",
				"entries for symbol", names(brackets[i]), "-", brackets[i]))
		}
	}

	output <- output[order(output[, 1]), ]
	col <- output[, 5]
	output <- output[, -5]
	output <- as.helix(output, length(chars))
	rownames(output) <- NULL

	if (!all(is.na(palette))) {
		palette <- rep(palette, length.out = 8)
		output$col <- palette[col]
	}

	return(output);
}


#' @export
#'
# Takes in a helix data structure and returns the corresponding structure
# Will allow up to 8 levels of pseudo-knots, and will REMOVE all conflicting
# basepairs, both done in a greedy fashion.
helixToVienna <- function(helix) {
	if(!is.helix(helix)) {
		stop("Invalid input detected, aborting")
	}
	if (any(helix$length > 1)) {
		warning("Expanding helix to basepairs")
		helix <- expandHelix(helix)
	}
	conflict <- isConflictingHelix(helix)
	if (any(conflict)) {
		warning("Removing conflicting basepairs")
		helix <- helix[!conflict, ]
	}
	if (nrow(helix) == 0) {
		stop("No basepairs in helix")
	}

	lbrac <- c("(", "<", "{", "[", "A", "B", "C", "D")
	rbrac <- c(")", ">", "}", "]", "a", "b", "c", "d")

	group <- unknottedGroups(helix)
	if (max(group) > length(lbrac)) {
		warning("Too much pseudoknots! Only returning first 8 levels")
	}
	str <- rep(".", attr(helix, "length"))
	for(i in 1:min(length(lbrac), max(group))) {
		str[helix$i[which(group == i)]] <- lbrac[i]
		str[helix$j[which(group == i)]] <- rbrac[i]
	}
	return(paste(str, collapse = ""))
}


#' @export
#'
# TODO: test with conflicting basepairs
helixToBpseq <- function(helix) {
	if (!is.helix(helix)) {
		stop("Not a valid helix data.frame, aborting")
	}
	if (is.null(attr(helix, "sequence"))) {
		stop("Missing nucleotide sequence, aborting")
	}
	output <- data.frame(i = as.integer(1:attr(helix, "length")),
		base = unlist(strsplit(attr(helix, "sequence"), "")))
	expand <- expandHelix(helix)[, c("i", "j")]
	output <- merge(output, expand, all = TRUE)
	output$j[is.na(output$j)] <- 0
	return(output)
}


#' @export
#'
helixToConnect <- function(helix) {
	bpseq <- helixToBpseq(helix)
	output <- data.frame(i = bpseq$i, base = bpseq$base, left = bpseq$i - 1,
		right = bpseq$i + 1, j = bpseq$j, same = bpseq$i)
	return(output)
}
