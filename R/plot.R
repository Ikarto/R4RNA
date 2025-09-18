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
blankPlot <- function(width, top, bottom, pad = c(0, 0, 0, 0),
		scale = TRUE, scale.lwd = 1, scale.col = "#DDDDDD", scale.cex = 1,
		debug = FALSE, png = NA, pdf = NA, factor = ifelse(!is.na(png), 8, 1/9),
		no.par = FALSE, ...) {
	x <- width + pad[2] + pad[4]
	# x <- max(x, 115) # width of longest legend
	y <- top - bottom + pad[1] + pad[3]
	if (!is.na(png)) {
		png(png, x * factor, y * factor)
	} else {
		if (!is.na(pdf)) {
			pdf(pdf, x * factor, y * factor)
		}
	}
	if (!no.par) { par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), ...) }
	plot(c(0 - pad[2], width + pad[4]), c(bottom - pad[1], top + pad[3]), type = "n", axes = F,
		xlab = "", ylab = "", asp = 1)
	if (scale) {
		plotScale(width, top, bottom, col = scale.col, lwd = scale.lwd, cex = scale.cex, xpd = TRUE)
	}
	if (debug) {
		print(paste("width", width, sep = ": "))
		print(paste("top", top, sep = ": "))
		print(paste("bottom", bottom, sep = ": "))
		print(paste("left", -pad[2], sep = ": "))
		print(paste("right", pad[4], sep = ": "))
		abline(h = c(top, bottom), col = "red", lty = 2)
		abline(v = c(-pad[2], width + pad[4]), col = "red")
		abline(h = 0, col = "green")
		abline(v = width, col = "green")
	}
}

plotScale <- function(width, top, bottom, ...) {
	pt <- strheight("0")
	bars <- pretty(c(1, width), n = 10)
	bars <- bars[which(bars <= width)]
	segments(bars, bottom - 1, y1 = top + 1.5 * pt, ...)
	text(bars, top + 1.5 * pt, bars, adj = c(-0.25, 1), ...)
}


#' @export
#'
maxHeight <- function(helix) {
	if(nrow(helix) == 0) {
		return(0)
	} else {
		return(max(abs(helix$i - helix$j)) / 2)
	}
}


arc <- function (c, r, v, theta, ...) {
    angles <- anglesArc(v, theta)
    seqang <- seq(angles[1], angles[2], length = 100)
    x <- c[1] + r * cos(seqang)
    y <- c[2] + r * sin(seqang)
    lines(x, y, ...)
}

anglesArc <- function (v, theta) {
    theta.OX <- ifelse(v[2] >= 0, acos(v[1]), 2 * pi - acos(v[1]))
    angs <- c(theta.OX - theta, theta.OX + theta)
    return(angs)
}

plotArrow <- function(start, y = 0, length = start * 0.03) {
	# NOTE: R has an automatic 4% margin for pretty axes, length is 3% default
	h <- (length / 4)
	polygon(c(start, start + length, start), c(h, 0, -h) + y, col = 1)
}

priority <- c("GC" = 1, "CG" = 2, "UA" = 3, "AU" = 4, "GU" = 5, "UG" = 6)


getCovarianceColours <- function(msa, helix) {
  cols <- matrix(5, nrow = (length(msa)), ncol = nchar(msa[1]))
  bases <- strsplit(msa, "")
  chars <- matrix(nrow = length(bases), ncol = length(bases[[1]]))
  if (nrow(helix) > 0 & is.helix(helix)) {
    for (i in 1:length(msa)) {
      bases[[i]] <- toupper(bases[[i]])
      bases[[i]] <- sub("T", "U", bases[[i]])
      chars[i, ] <- bases[[i]]
      col <- getSequenceColour(msa[i], helix)
      cols[i, ] <- col
    }
    if (length(msa) == 1) {
      cols[1, which(cols[1, ] == 2)] <- 1
    } else {
      for (i in 1:nrow(helix)) {
        # Colours conservation relative to the most common valid basepair
        pos <- as.integer(helix[i, c("i", "j")])
        bases <- chars[, pos]
        pairs <- apply(bases, 1, paste, collapse = "")
        counts <- table(pairs[which(cols[, pos[1]] == 2)])
        if(length(counts) > 0) {
          max <- names(priority)[min(priority[names(which(counts ==
                                                            max(counts)))])]
          cols[which(pairs == max), pos] <- 1

          # One-sided conservation colouring
          maxchar <- unlist(strsplit(max, ""))
          blues <- which(cols[, pos[1]] == 2)
          hits <- blues[which(bases[blues, 1] == maxchar[1])]
          cols[hits, pos[1]] <- 1
          cols[hits, pos[2]] <- 3
          hits <- blues[which(bases[blues, 2] == maxchar[2])]
          cols[hits, pos[2]] <- 1
          cols[hits, pos[1]] <- 3
        }
      }
    }
  }
  return(cols)
}

getCovarianceColoursMltp <- function(msa, helix) {
	cols <- matrix(5, nrow = (length(msa)), ncol = Biostrings::nchar(msa[1]))
	bases <- as.matrix(msa)
	chars <- matrix(nrow = nrow(bases), ncol = ncol(bases))
	if (nrow(helix) > 0 & is.helix(helix)) {
		for (i in 1:length(msa)) {
			bases[i,] <- toupper(bases[i,])
			bases[i,] <- sub("T", "U", bases[i,])
			chars[i, ] <- bases[i,]
			col <- getSequenceColour(msa[i], helix)
			cols[i, ] <- col
		}
		if (length(msa) == 1) {
			cols[1, which(cols[1, ] == 2)] <- 1
		} else {
			for (i in 1:nrow(helix)) {
				# Colours conservation relative to the most common valid basepair
				pos <- as.integer(helix[i, c("i", "j")])
				bases <- chars[, pos]
				pairs <- apply(bases, 1, paste, collapse = "")
				counts <- table(pairs[which(cols[, pos[1]] == 2)])
				if(length(counts) > 0) {
					max <- names(priority)[min(priority[names(which(counts ==
						max(counts)))])]
					cols[which(pairs == max), pos] <- 1

					# One-sided conservation colouring
					maxchar <- unlist(strsplit(max, ""))
					blues <- which(cols[, pos[1]] == 2)
					hits <- blues[which(bases[blues, 1] == maxchar[1])]
					cols[hits, pos[1]] <- 1
					cols[hits, pos[2]] <- 3
					hits <- blues[which(bases[blues, 2] == maxchar[2])]
					cols[hits, pos[2]] <- 1
					cols[hits, pos[1]] <- 3
				}
			}
		}
	}
	return(cols)
}

getBaseColours <- function(msa) {
	msa <- toupper(msa)
	msa <- gsub("[.]", "-", msa)
	msa <- gsub("T", "U", msa)
	msa <- gsub("[^ACGU-]", "N", msa)
	cols <- matrix(7, nrow = length(msa), ncol = Biostrings::nchar(msa[1]))
	map <- c("A" = 1, "U" = 2, "G" = 3, "C" = 4, "-" = 5, "N" = 6)
	for(i in 1:length(msa)) {
		s <- msa[i]
		for(j in 1:nchar(s)) {
			col <- map[as.character(toupper(substr(s, j, j)))]
			if(!is.na(col)) { cols[i, j] <- col }
		}
	}
	return(cols)
}

getSequenceColour <- function(msa, helix) {
	msa <- toupper(msa)
	msa <- gsub("[.]", "-", msa)
	msa <- gsub("T", "U", msa)
	msa <- gsub("[^ACGU-]", "N", msa)
	cols <- rep(5, nchar(msa))
	chars <- unlist(strsplit(msa, ""))
	cols[which(chars == '-')] <- 6
	cols[which(chars == 'N')] <- 7
	for(i in 1:nrow(helix)) {
		pos <- as.integer(helix[i, c("i", "j")])
		col <- bpCols[chars[pos[1]], chars[pos[2]]]
		cols[pos[1]] <- col[[1]][1]
		cols[pos[2]] <- col[[1]][2]
	}
	return(cols)
}

bpCols <- matrix(list(c(4, 4)), ncol = 6, nrow = 6)
rownames(bpCols) <- c("G", "C", "A", "U", "-", "N")
colnames(bpCols) <- c("G", "C", "A", "U", "-", "N")
bpCols["A", "U"] <- list(c(2, 2))
bpCols["U", "A"] <- list(c(2, 2))
bpCols["G", "C"] <- list(c(2, 2))
bpCols["C", "G"] <- list(c(2, 2))
bpCols["G", "U"] <- list(c(2, 2))
bpCols["U", "G"] <- list(c(2, 2))
bpCols["-", "G"] <- list(c(6, 4))
bpCols["-", "C"] <- list(c(6, 4))
bpCols["-", "A"] <- list(c(6, 4))
bpCols["-", "U"] <- list(c(6, 4))
bpCols["G", "-"] <- list(c(4, 6))
bpCols["C", "-"] <- list(c(4, 6))
bpCols["A", "-"] <- list(c(4, 6))
bpCols["U", "-"] <- list(c(4, 6))
bpCols["-", "-"] <- list(c(6, 6))
bpCols["N",    ] <- list(c(7, 7))
bpCols[   , "N"] <- list(c(7, 7))
bpCols["N", "-"] <- list(c(7, 6))
bpCols["-", "N"] <- list(c(6, 7))
#################################
#bpCols["A", "T"] <- list(c(2, 2))
#bpCols["T", "A"] <- list(c(2, 2))
#bpCols["G", "T"] <- list(c(2, 2))
#bpCols["T", "G"] <- list(c(2, 2))
#bpCols["-", "T"] <- list(c(6, 4))
#bpCols["T", "-"] <- list(c(4, 6))

#' @export
#'
plotArc <- function (i, j, y = 0, flip = FALSE, shape = "circle", ...) {
  i <- as.numeric(i)
  j <- as.numeric(j)
  if (shape == "circle") {
    center <- c(mean(c(i, j)), y)
    radius <- abs(diff(c(i, j)))/2
    vector <- c(0, ifelse(flip, -1, 1))
    theta <- pi/2
    arc(center, radius, vector, theta, ...)
  }else {
    width <- abs(i - j)
    height <- (width/2) + y
    height <- ifelse(flip, -height, height)
    if (shape == "triangle") {
      center <- (i + j)/2
      lines(c(i, center, j), c(y, height, y), ...)
    }else {
      if (shape == "square") {
        lines(c(i, i, j, j), c(y, height, height, y),...)
      }else{
        if (shape == "heatmap") {
          center <- (i + j)/2
          polygon(c(center-0.5,center,center+0.5,center),
                  c(height,height-0.5,height,height+0.5)
                  ,border = NA,...)
        }
      }
    }
  }
}



#' @export
#'
plotArcs <- function(i, j, length, y = 0, flip = FALSE, shape = "circle", ...) {
	i <- as.numeric(i)
	j <- as.numeric(j)
	for (n in 0:(length - 1)) {
		plotArc(i + n, j - n, y, flip, shape = shape, ...)
	}
}


#' @export
#'
plotHelix <- function(helix, y = 0, flip = FALSE, line = FALSE, arrow = FALSE,
		add = FALSE, shape = "circle", ...) {
	if (!is.helix(helix)) {
		stop("Input not a valid helix data frame, aborting")
	}
	if (nrow(helix) > 0) {
		if (is.null(helix$flip)) {
			helix$flip <- flip
		}
		if (is.null(helix$y)) {
			helix$y <- y
		}
		if (is.null(helix$shape)) {
			helix$shape <- shape
		}
		style <- NULL
		width <- attr(helix, "length")
		if(ncol(helix) > 4) {
			types <- lapply(helix, typeof)
			names <- names(helix)
			for (i in 5:ncol(helix)) {
				helix[names[i]] <- paste(names[i], " = ", "as.", types[i], "('",
					unlist(helix[i]), "')", sep = "")
			}
			helix$cmd <- apply(data.frame(helix[, 5:ncol(helix)]), 1, paste, collapse = ", ")
			helix <- helix[, c("i", "j", "length", "value", "cmd")]
		}
		if (!add) {
			blankPlot(width, maxHeight(helix), y, ...)
		}

		# Reversed so top helices appear on top of lower helices
		for (n in nrow(helix):1) {
			if (!is.null(helix$cmd)) {
				cmd <- paste("plotArcs(helix$i[n], helix$j[n], helix$length[n],", helix$cmd[n], ")")
				eval(parse(text = cmd))
			} else {
				plotArcs(helix$i[n], helix$j[n], helix$length[n], y = y, flip = flip,
					shape = shape)
			}

		}
		if (line) {
			lines(c(0.5, width + 0.5), c(y, y), col = 1)
		}

		if (arrow) {
			plotArrow(width + 0.5, y)
		}
	} else {
		warning("Input structure is empty, no helices plotted")
	}

}



#' @export
#'
plotDoubleHelix <- function(top, bot, line = TRUE, arrow = FALSE, add = FALSE,
		...) {
	if (!is.helix(top) || !is.helix(bot)) {
		stop("Input not a valid helix data frame, aborting")
	}
	width <- max(attr(top, "length"), attr(bot, "length"))
	if (!add) {
		blankPlot(width, maxHeight(top), -maxHeight(bot), ...)
	}
	plotHelix(top, add = TRUE, ...)
	plotHelix(bot, flip = TRUE, add = TRUE, ...)

	if (line) {
		lines(c(0.5, width + 0.5), c(0, 0), col = 1)
	}

	if (arrow) {
		plotArrow(width + 0.5)
	}
}



#' @export
#'
plotOverlapHelix <- function(predict, known, miss = "black", line = TRUE,
	arrow = FALSE, add = FALSE, any = FALSE, ...) {
	if (!is.helix(predict) || !is.helix(known)) {
		stop("Input not a valid helix data frame, aborting")
	}
	width <- max(attr(predict, "length"), attr(known, "length"))

	predict$flip <- FALSE
	overlap <- isOverlappingHelix(predict, known, any = any)
	predict$flip[which(!overlap)] <- TRUE

	misses <- known[which(!isOverlappingHelix(known, predict)), ]

	if (!add) {
		top <- max(maxHeight(predict[which(predict$flip == FALSE), ]),
			maxHeight(misses))
		bot <- maxHeight(predict[which(predict$flip == TRUE), ])
		blankPlot(width, top, -bot, ...)
	}

	if (nrow(misses) > 0) {
		plotHelix(misses, col = miss, add = TRUE, ...)
	}
	plotHelix(predict, add = TRUE, ...)

	if (line) {
		lines(c(0.5, width + 0.5), c(0, 0), col = 1)
	}
	if (arrow) {
		plotArrow(width + 0.5)
	}
}



plotCovarianceLine <- function(x, y, cols, ...) {
	x <- x - 0.5
	start <- x:(x + length(cols) - 1)
	end <- start + 1
	segments(start, y, end, col = cols, lend = 1, ...)
}

plotCovarianceGrid <- function(x, y, cols, height = 1.25, ...) {
	x <- x - 0.5
	start <- x:(x + length(cols) - 1)
	end <- start + 1
	top <- y
	bot <- y + height
	rect(start, bot, end, top, col = cols, ...)
}

plotCovarianceText <- function(x, y, text, height = 1.25, ...) {
	# font = 2, family = "sans"
	chars <- unlist(strsplit(text, ""))
	start <- 1:length(chars) + x - 1
	# 'font': c(san-serif, bold, italics, bold italics, serif)
	# 'family': "mono", "serif", "sans"
	text(start, y + height / 2, chars, ...)
}

plotCovarianceSpecies <- function(x, y, text, height = 1.25, ...) {
 	# font = 3, family = "mono", adj = 0
	text(-x, y + height / 2, text, adj = 0, ...)
}



#' @export
#'
plotOverlapCovariance <- function(predict.helix, known.helix, msa, bot.msa = TRUE,
		any = FALSE, miss = "black", add = FALSE, grid = FALSE, species = 0,
		legend = TRUE, pad = c(0, 0, 0, 0), ...) {
	if (!is.helix(predict.helix) || !is.helix(known.helix)) {
		stop("One of input not a valid helix data frame, aborting")
	}
	width <- max(nchar(msa))
	height <- ifelse(grid, 1.25, 1)

	overlap <- isOverlappingHelix(predict.helix, known.helix, any = any)
	tp <- predict.helix[overlap, ]
	fp <- predict.helix[!overlap, ]
	fn <- known.helix[!isOverlappingHelix(known.helix, predict.helix), ]
	if (nrow(fn) > 0) {
		fn$col <- miss
	}

	offset <- (length(msa) + 0.5) * height
	up <- max(maxHeight(tp), maxHeight(fn)) + offset
	down <- maxHeight(fp) + ifelse(bot.msa, offset, 0)

	if (!add) {
		if (species > 0) {
			pad[2] <- pad[2] + species
		}
		suppressWarnings(blankPlot(width, up, -down, pad, ...))
	}

	if (bot.msa) {
		plotCovariance(msa, fp, add = TRUE, flip = TRUE, grid = grid, species = species, legend = FALSE, y = -height/2, ...)
	} else {
		plotHelix(fp, add = TRUE, flip = TRUE, y = height/2, ...)
	}
	if (nrow(fn) > 0) {
		plotHelix(fn, add = TRUE, y = offset - ifelse(grid, 0, 1), ...)
	}
	plotHelix(tp, add = TRUE, y = offset - ifelse(grid, 0, 1), ...)
	plotCovariance(msa, known.helix, add = TRUE, arc = FALSE, grid = grid, species = species, legend = legend, y = height/2, ...)
}




#' @export
#'
plotDoubleCovariance <- function(top.helix, bot.helix, top.msa, bot.msa =
		top.msa, add = FALSE, grid = FALSE, species = 0, legend = TRUE,
		pad = c(0, 0, 0, 0), ...) {
	if (!is.helix(top.helix) || !is.helix(bot.helix)) {
		stop("One of input not a valid helix data frame, aborting")
	}
	width <- max(nchar(c(top.msa, bot.msa)), na.rm = TRUE)
	height <- ifelse(grid, 1.25, 1)

	up <- maxHeight(top.helix) + (length(top.msa) + 0.5) * height
	down <- maxHeight(bot.helix) + (length(bot.msa) + 0.5) * height

	if (!add) {
		if (species > 0) {
			pad[2] <- pad[2] + species
		}
		suppressWarnings(blankPlot(width, up, -down, pad, ...))
	}

	if (any(is.na(bot.msa))) {
		plotHelix(bot.helix, flip = TRUE, add = TRUE, y = height/2, ...)
	} else {
		plotCovariance(bot.msa, bot.helix, flip = TRUE, add = TRUE, grid = grid, species = species, legend = FALSE, y = -height/2, ...)
	}
	plotCovariance(top.msa, top.helix, add = TRUE, grid = grid, species = species, legend = legend, y = height/2, ...)
}



#' @export
#'
plotCovariance <- function(msa, helix, arcs = TRUE, add = FALSE, grid = FALSE,
		text = FALSE, legend = TRUE, species = 0, base.colour = FALSE, palette =
		NA, flip = FALSE, grid.col = "white", grid.lwd = 0, text.cex = 0.5,
		text.col = "white", text.font = 2, text.family = "sans", species.cex =
		0.5, species.col = "black", species.font = 2, species.family = "mono",
		shape = "circle", y = 0, conflict.lty = 2, conflict.col = NA,
		pad = c(0, 0, 0, 0), ...) {
	if (!is.helix(helix)) {
		stop("Input not a valid helix data frame, aborting")
	}
	x <- 1
	width <- max(nchar(msa))
	height <- ifelse(grid, 1.25, 1)

	if (any(nchar(msa) != width)) {
		warning("WARNING: Padding sequences of unequal length with gaps")
		unequal <- which(nchar(msa) != width)
		for (i in 1:length(unequal)) {
			msa[unequal[i]] <- paste(msa[unequal[i]], paste(rep("-",
				width - nchar(msa[unequal[i]])), collapse = ""), sep = "")
		}
	}

	if (all(is.na(palette))) {
		if (base.colour) {
			palette <- c("#E41A1C", "#4DAF4A", "#FF7F00", "#377EB8", "#BDBDBD",
				"#984EA3")
		} else {
			palette <- c("#00A651", "#0072BC", "#00B9F2", "#F15A22", "#231F20",
				"#AAAAAA", "#DA6FAB")
		}
	} else {
		if (length(palette) < ifelse(base.colour, 5, 7)) {
			palette <- rep(palette, length.out = ifelse(base.colour, 5, 7))
		}
	}

	if (flip) {
		up <- y
		down <- y - maxHeight(helix) - length(msa) * height
	} else {
		down <- y
		y <- y + length(msa) * height
		up <- maxHeight(helix) + y - ifelse(grid, 0, 1)
	}

	if (!add) {
		if (species > 0) {
			pad[2] <- pad[2] + species
		}
		blankPlot(width, up, down, pad, ...)
	}

	conflict <- NULL
	if (!is.na(conflict.col) || !is.na(conflict.lty)) {
		conflicting <- isConflictingHelix(helix)
		conflict <- helix[conflicting, ]
		helix <- helix[!conflicting, ]

		if (nrow(conflict) > 0) {
			if (!is.na(conflict.col)) {
				conflict$col <- conflict.col
			}
			if (!is.na(conflict.lty)) {
				conflict$lty <- conflict.lty
			}
		} else {
			conflict <- NULL
		}
	}

	if (arcs & nrow(helix) != 0) {
		if (flip) {
			if (!is.null(conflict)) {
				plotHelix(conflict, y = y - (length(msa)) * height, add = TRUE,
					flip = flip, ...)
			}
			plotHelix(helix, y = y - (length(msa)) * height, add = TRUE,
				flip = flip, ...)
		} else {
			if (!is.null(conflict)) {
				plotHelix(conflict, y = y - ifelse(grid, 0, 1), add = TRUE, ...)
			}
			plotHelix(helix, y = y - ifelse(grid, 0, 1), add = TRUE, ...)
		}
	}

	if (base.colour) {
		cols <- getBaseColours(msa)
	} else {
		cols <- getCovarianceColours(msa, expandHelix(helix))
	}

	used <- sort(unique(c(cols)))
	dim <- dim(cols)
	cols <- palette[cols]
	dim(cols) <- dim

	for(i in 1:length(msa)) {
		if (grid) {
			plotCovarianceGrid(x, y - height * i, cols[i, ], height = height,
				border = grid.col, lwd = grid.lwd)
			if (text) {
				plotCovarianceText(x, y - height * i, msa[i], height = height,
					cex = text.cex, font = text.font, family = text.family,
					col = text.col)
			}
			if (species > 0) {
				plotCovarianceSpecies(species - x, y - height * i,
					names(msa[i]), height = height, cex = species.cex,
					font = species.font, family = species.family,
					col = species.col)
			}
		} else {
			plotCovarianceLine(x, y - i, cols[i, ])
		}
	}

	if (legend) {
		if (base.colour) {
			text <- c("A", "U", "G", "C", "-", "?")
		} else {
			text <- c("Conservation", "Covariation", "One-sided", "Invalid",
				"Unpaired", "Gap", "Ambiguous")
		}
		legend("bottom", text[used], border = NA, fill = palette[used],
			horiz = TRUE, bty = "n", xpd = NA)
	}
}



#' @export
#'
colourByUnknottedGroups <- function(helix, cols, get = TRUE) {
	if (!is.helix(helix)) {
		stop("Not a valid helix data frame, aborting")
	}
	if (missing(cols)) {
		cols <- defaultPalette()
	}
	expanded <- expandHelix(helix)
	group <- as.factor(unknottedGroups(expanded))
	if (length(levels(group)) > length(cols)) {
		warning(paste("Number of groups (", length(levels(group)), ") greater than colours (", length(cols) ,"), some groups will be colourless", sep = ""))
	}
	output <- cols[as.integer(group)]
	legend <- table(output)[cols]
	legend[is.na(legend)] <- 0
	if (get) {
		expanded$col <- output
		attr(expanded, "legend") <- paste(legend, "/", sum(legend), sep = "")
		attr(expanded, "fill") <- cols
		return(expanded)
	} else {
		attr(output, "legend") <- paste(legend, "/", sum(legend), sep = "")
		attr(output, "fill") <- cols
		return(output)
	}
}



#' @export
#'
colourByBasepairFrequency <- function(helix, cols, get = TRUE) {
	if (!is.helix(helix)) {
		stop("Not a valid helix data frame, aborting")
	}
	if (missing(cols)) {
		cols <- rev(defaultPalette())
	}
	freq <- colourByValue(basepairFrequency(helix), cols, get = TRUE)
	if (get) {
		return(freq)
	} else {
		output <- freq$col
		attr(output, "legend") <- attr(freq, "legend")
		attr(output, "fill") <- attr(output, "fill")
		return(output)
	}
}



#' @export
#'
colourByValue <- function(helix, cols, breaks, get = FALSE, log = FALSE, include.lowest = TRUE, ...) {
	if (!is.helix(helix)) {
		stop("Not a valid helix data frame, aborting")
	} else {
		if (all(is.na(helix$value))) {
			stop("Cowardly refusing to deal with NA/NaN values")
		}
	}
	if (missing(cols)) {
		cols <- defaultPalette()
	}
	if (missing(breaks)) {
		breaks <- seq(min(helix$value, na.rm = TRUE),
			max(helix$value, na.rm = TRUE), length.out = length(cols) + 1)
	} else {
		if (length(breaks) == 1 && breaks == 1) {
			breaks <- c(min(helix$value, na.rm = TRUE),
				max(helix$value, na.rm = TRUE))
		}
	}
	if (log) {
		test <- helix$value[helix$value > 0]
		start <- logfloor(min(test, na.rm = TRUE))
		end <- logceiling(max(test, na.rm = TRUE))
		breaks <- c(min(helix$value, na.rm = TRUE), logseq(start, end)[-1])
	}
	levels <- cut(helix$value, breaks = breaks, labels = NULL, ordered_result = FALSE, include.lowest = include.lowest, ...)
	legend <- levels(levels)
	if (length(legend) > length(cols)) {
		warning(paste("Number of intervals (", length(legend), ") greater than colours (", length(cols) ,"), some intervals will be colourless", sep = ""))
	}
	fill <- cols[1:length(legend)]
	output <- cols[as.numeric(levels)]
	attr(output, "legend") <- legend
	attr(output, "fill") <- fill
	if (get) {
		helix$col <- output
		attr(helix, "legend") <- legend
		attr(helix, "fill") <- fill
		return(helix)
	} else {
		return(output)
	}
}


#' @export
#'
defaultPalette <- function() {
	return(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))
}



#' @export
#'
colourByCount <- function(helix, cols, counts, get = FALSE) {
	if (!is.helix(helix)) {
		stop("Not a valid helix data frame, aborting")
	}
	if (missing(cols)) {
		cols <- defaultPalette()
	}
	if (missing(counts)) {
		output <- as.character(cut(seq_len(nrow(helix)), breaks = seq(1,
			nrow(helix), length.out = length(cols) + 1), include.lowest = TRUE,
			labels = cols))
	} else {
		if (length(cols) != length(counts)) {
			stop(paste("Length of counts (", length(counts),
				") does not match length of cols (", length(cols), ")",
				sep = ""))
		}
		output <- as.character(rep(cols, counts))
	}
	legend <- table(output)[cols]
	legend[is.na(legend)] <- 0
	if (get) {
		helix$col <- NA
		n <- min(length(output), nrow(helix))
		helix$col[1:n] <- output[1:n]
		attr(helix, "legend") <- paste(legend, "/", sum(legend), sep = "")
		attr(helix, "fill") <- cols
		return(helix)
	} else {
		attr(output, "legend") <- paste(legend, "/", sum(legend), sep = "")
		attr(output, "fill") <- cols
		return(output)
	}
}

