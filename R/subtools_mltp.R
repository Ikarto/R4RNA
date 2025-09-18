## Copyright (C) 2019 Irmtraud M. Meyer,Volodymyr Tsybuslkyi, Daniel Lai, (www.e-rna.org)
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


###############################################################################
#getBaseColoursMod
###############################################################################
getBaseColoursMod <- function(msa) {
  msa <- toupper(msa)
  msa <- gsub("[.]", "-", msa)
  msa <- gsub("[^ACGUT-]", "N", msa)
  cols <- matrix(7, nrow = length(msa), ncol = Biostrings::nchar(msa[1]))
  map <- c("A" = 1, "U" = 2, "G" = 3, "C" = 4, "-" = 5, "N" = 6,"T" = 7)
  for(i in 1:length(msa)) {
    s <- msa[i]
    for(j in 1:nchar(s)) {
      col <- map[as.character(toupper(substr(s, j, j)))]
      if(!is.na(col)) { cols[i, j] <- col }
    }
  }
  return(cols)
}

###############################################################################
#update helix with coordinates
###############################################################################
update.helix.coord.double <- function(helix,uniq,start.uniq,end.uniq,mode = 1){
  if(mode ==1){
    for(k in start.uniq:end.uniq) {
      helix[ent.i %like% uniq[k,group],i.coord := i+uniq[k,start]]
      helix[ent.i %like% uniq[k,group],ent.i.y := uniq[k,y]]
      helix[ent.j %like% uniq[k,group],j.coord := j+uniq[k,start]]
      helix[ent.j %like% uniq[k,group],ent.j.y := uniq[k,y]]
    }
  }
  if(mode==2){
    for(k in start.uniq:end.uniq) {
      helix[ent.i %like% uniq[k,group],i.coord := uniq[k,end]-i+1]
      helix[ent.i %like% uniq[k,group],ent.i.y := uniq[k,y]]
      helix[ent.j %like% uniq[k,group],j.coord := uniq[k,end]-j+1]
      helix[ent.j %like% uniq[k,group],ent.j.y := uniq[k,y]]
    }
  }
  if(mode==3){
    for(k in 1:uniq[,.N]) {
      if(k==1){
        helix[ent.i %like% uniq[k,group],i.coord := uniq[k,start]+i]
        helix[ent.i %like% uniq[k,group],ent.i.y := uniq[k,y]]
        helix[ent.j %like% uniq[k,group],j.coord := uniq[k,start]+j]
        helix[ent.j %like% uniq[k,group],ent.j.y := uniq[k,y]]
      }else{
        helix[ent.i %like% uniq[k,group],i.coord := uniq[k,end]-i+1]
        helix[ent.i %like% uniq[k,group],ent.i.y := uniq[k,y]]
        helix[ent.j %like% uniq[k,group],j.coord := uniq[k,end]-j+1]
        helix[ent.j %like% uniq[k,group],ent.j.y := uniq[k,y]]
      }
    }
  }
  return(helix)
}



###############################################################################
#calculation for single line mod
###############################################################################
prep.uniq.single <- function(uniq,dist.part,full.len = 0){
  #set index N for location in plot
  index <- uniq[,.I]
  uniq[,N:=index]
  len <- uniq[,as.numeric(width)]
  uniq[,width:=len]
  #distance between entities
  dist.between <- as.integer(dist.part*mean(uniq[,width]))
  #define start and end pos
  if(full.len == 0){
    full.len <- sum(uniq[,width]) + (uniq[,.N]-1)*dist.between
  }
  for (i in 1:uniq[,.N]) {
    uniq[i,start := full.len-(uniq[,.N]-i)*dist.between-
           sum(uniq$width[i:uniq[,.N]])]
  }
  uniq[,end := start+width]
  attr(uniq,"full.length") <- full.len
  attr(uniq,"dist.between") <- dist.between
  return(uniq)
}


###############################################################################
#plot bases for single mode
###############################################################################
plot.bases.single <- function(uniq,col.line,y=0){
  if("y" %in% colnames(uniq)){
    for(i in 1:nrow(uniq)){
      lines(c(uniq[i,start],uniq[i,end]),c(uniq[i,y],uniq[i,y]),col = col.line)
    }
  }else{
    for(i in 1:nrow(uniq)){
      lines(c(uniq[i,start],uniq[i,end]),c(y,y),col = col.line)
    }
  }
}

###############################################################################
#plot arrow in the end of bases
###############################################################################
plot.arrows.single <- function(uniq,dist.between,col.arrow = "black",mode = 1,y=0){
  if(mode ==1){
    if("y" %in% colnames(uniq)){
      for (i in 1:uniq[,.N]) {
        x.coord <- c(uniq[i,end],uniq[i,end],uniq[i,end+dist.between/8])
        y.coord <- c(uniq[i,y]-1,uniq[i,y]+1,uniq[i,y])
        polygon(x.coord, y.coord,col = col.arrow)
      }
    }else{
      for (i in 1:uniq[,.N]) {
        x.coord <- c(uniq[i,end],uniq[i,end],uniq[i,end]+dist.between/8)
        y.coord <- c(y-1,y+1,y)
        polygon(x.coord, y.coord,col = col.arrow)
      }
    }
  }
  if(mode ==2){
    if("y" %in% colnames(uniq)){
      for (i in 1:uniq[,.N]) {
        if(i==1){
          x.coord <- c(uniq[i,end],uniq[i,end],uniq[i,end+dist.between/8])
          y.coord <- c(uniq[i,y]-1,uniq[i,y]+1,uniq[i,y])
          polygon(x.coord, y.coord,col = col.arrow)
        }else{
          x.coord <- c(uniq[i,start],uniq[i,start],uniq[i,start-dist.between/8])
          y.coord <- c(uniq[i,y]-1,uniq[i,y]+1,uniq[i,y])
          polygon(x.coord, y.coord,col = col.arrow)
        }
      }
    }else{
      if(i==1){
        x.coord <- c(uniq[i,end],uniq[i,end],uniq[i,end+dist.between/8])
        y.coord <- c(y-1,y+1,y)
        polygon(x.coord, y.coord,col = col.arrow)
      }else{
        x.coord <- c(uniq[i,start],uniq[i,start],uniq[i,start-dist.between/8])
        y.coord <- c(y-1,y+1,y)
        polygon(x.coord, y.coord,col = col.arrow)
      }
    }
  }

}


###############################################################################
#update.uniq.by.group
###############################################################################
update.uniq.by.group.single <- function(uniq,top.name,sort){
  #check if name and sorting is specified
  options(warn = -1)

  if(top.name != "FALSE"){
    if(uniq[1,group]!=top.name){
      if(sort){
        uniq <- rbind(uniq[group %like% top.name],
                      uniq[!(group %like% top.name)][order(group)])
      }else{
        uniq <- rbind(uniq[group %like% top.name],
                      uniq[!(group %like% top.name)])
      }
    }
  }
  if(sort & top.name == "FALSE"){
    uniq <- uniq[order(group)]
  }
  return(uniq)
}



###############################################################################
#update.helix.coord
###############################################################################
update.helix.coord.single <- function(helix,uniq){
  #reduce data.table warning
  options(warn = -1)

  for(k in 1:uniq[,.N]) {
    helix[ent.i %like% uniq[k,group],i.coord := i+uniq[k,start]]
    helix[ent.j %like% uniq[k,group],j.coord := j+uniq[k,start]]
  }
  return(helix)
  options(warn = 1)
}

###############################################################################
#max.height.single
###############################################################################
max.height.single <- function(helix,full.len,stable) {
  if(stable){
    max.height <- full.len/2
  }else{
    max.height <- max(helix[,abs(j.coord-i.coord)])/2
  }
}


###############################################################################
#plot lines as trans-interactions
###############################################################################
plot.lines.single <- function(helix,col="black",lty = 1,lwd=1){
  if("col" %in% colnames(helix)){
    for (k in 1:helix[,.N]) {
      lines(x = c(helix[k,i.coord],helix[k,j.coord]), y = c(helix[k,ent.i.y],helix[k,ent.j.y]),col=helix[k,col],lty=lty,lwd=lwd)
    }
  }else{
    for (k in 1:helix[,.N]) {
      lines(x = c(helix[k,i.coord],helix[k,j.coord]), y = c(helix[k,ent.i.y],helix[k,ent.j.y]),col=col,lty=lty,lwd=lwd)
    }
  }
}



###############################################################################
#define max height for helix
###############################################################################
def.max.height.helix <- function(helix,top.name = "FALSE",
                                 sort = TRUE,dist.part = 0.2,
                                 stable = FALSE){
  helix <-  copy(helix)
  uniq <- attr(helix,"length")
  uniq <- update.uniq.by.group.single(uniq,top.name,sort)
  uniq <- prep.uniq.single(uniq,dist.part)
  full.len <- attr(uniq,"full.length")
  dist.between <- attr(uniq,"dist.between")
  helix.processed <- update.helix.coord.single(helix,uniq)
  max.height <- max.height.single(helix.processed,full.len,stable)
  attr(max.height,"dist.between") <- dist.between
  attr(max.height,"uniq") <- uniq
  attr(max.height,"full.length") <-full.len
  return(max.height)
}


###############################################################################
#process readed fasta by readBStringSet() function
###############################################################################
bstring.to.data.table <- function(bstring,comp.seq = TRUE){
  if(comp.seq){
    #extract names and define entity group and sequence name identification
    #delimiter is always "."
    names <- as.vector(names(bstring))
    tmp <- regmatches(names, regexpr("[.]", names), invert = TRUE)
    group <- c()
    group <- sapply(tmp, function(x) group <- c(group,x[1]),USE.NAMES = FALSE)
    if(NA %in% group){
      stop("\nERROR: Sequence without entity name introduced")
    }
    tmp.seq.name <- c()
    tmp.seq.name <- sapply(tmp, function(x) tmp.seq.name <- c(tmp.seq.name,x[2]),USE.NAMES = FALSE)
    if(NA %in% tmp.seq.name){
      stop("\nERROR: Sequence without sequence name introduced")
    }
    tmp.seq.name <- regmatches(tmp.seq.name, regexpr(":", tmp.seq.name), invert = TRUE)
    seq.name <- c()
    seq.name <- sapply(tmp.seq.name, function(x) seq.name <- c(seq.name,x[1]),USE.NAMES = FALSE)
    rest <- c()
    rest <- sapply(tmp.seq.name, function(x) rest <- c(rest,x[2]),USE.NAMES = FALSE)
    if(NA %in% rest){
      message("\nSequences without rest part introduced.\nIt will not influence work of plotting\n")
    }
    rest <- as.data.table(rest)
    rest[which(is.na(rest)),] <- " "
    sequence <- as.character(bstring)
    width <- base::nchar(sequence)
    #create data.table
    fasta <- data.table(as.data.table(group), as.data.table(seq.name),as.data.table(rest),
                        as.data.table(sequence),as.data.table(width))
    return(fasta)
  }else{
    group <- c()
    name <- c()
    names <- as.vector(names(bstring))
    tmp <- regmatches(names, regexpr("[.]", names), invert = TRUE)
    group <- sapply(tmp, function(x) group <- c(group,x[1]),USE.NAMES = FALSE)
    name <- sapply(tmp, function(x) name <- c(name,x[2]),USE.NAMES = FALSE)
    if(NA %in% group){
      stop("\nERROR: Sequence without entity name introduced")
    }
    if(NA %in% name){
      message("\nSequences without name part introduced.\nIt will not influence work of plotting\n")
    }
    sequence <- as.character(bstring)
    width <- base::nchar(sequence)
    fasta <- data.table(as.data.table(group), as.data.table(name),
                        as.data.table(sequence),as.data.table(width))
    return(fasta)
  }
}


###############################################################################
#anglesArc
###############################################################################
anglesArc <- function (v, theta) {
  theta.OX <- ifelse(v[2] >= 0, acos(v[1]), 2 * pi - acos(v[1]))
  angs <- c(theta.OX - theta, theta.OX + theta)
  return(angs)
}

###############################################################################
#arc
###############################################################################
arcDouble <- function(c, r, v, theta,y1,y2, ...) {
  angles <- anglesArc(v, theta)
  seqang <- seq(angles[1], angles[2], length = 100)
  x <- c[1] + r * cos(seqang)
  y <- c[2] + r * sin(seqang)
  if(y1<y2){
    difference <- abs(diff(c(y1,y2)))/2
    y <- y + seq(difference,-difference,length.out = 100)
  }
  if(y1>y2){
    difference <- abs(diff(c(y1,y2)))/2
    y <- y + seq(-difference,difference,length.out = 100)
  }
  lines(x, y, ...)
}

###############################################################################
#plot lines as trans-interactions
###############################################################################
plotArcDouble <- function(i, j, x = 0, y1 = 0,y2 = 0, flip = FALSE, shape = "circle", ...){
  i <- as.numeric(i) + x
  j <- as.numeric(j) + x

  if (shape == "circle") {
    center <- c(mean(c(i, j)), mean(c(y1, y2)))
    radius <- abs(diff(c(i, j))) / 2
    vector <- c(0, ifelse(flip, -1, 1))
    theta <- pi/2
    arcDouble(center, radius, vector, theta,y1,y2, ...)
  }
  else {
    width <- abs(i - j)
    height <- (width / 2) + mean(c(y1, y2))
    height <- ifelse(flip, -height, height)
    if (shape == "triangle") {
      center <- (i + j)/ 2
      lines(c(i, center, j), c(y1, height, y2), ...)
    } else {
      if (shape == "square") {
        lines(c(i, i, j, j), c(y1, height, height, y2), ...)
      }
    }
  }
}

###############################################################################
#plot arcs for single mod
###############################################################################
plot.arcs.single <- function(helix,col = "black",flip = FALSE,y = 0,lty=1,shape="circle",lwd = 1){
  if("col" %in% colnames(helix)){
    if("ent.i.y" %in% colnames(helix)){
      for (k in 1:helix[,.N]) {
        if(helix[k,ent.i.y]!=helix[k,ent.i.y]){
          plotArcDouble(i = helix[k,i.coord],j =helix[k,j.coord] ,y1 = helix[k,ent.i.y],y2 = helix[k,ent.j.y],flip = flip,col = helix[k,col],lty =lty,shape=shape,lwd = lwd)
        }else{
          plotArc(i = helix[k,i.coord],j =helix[k,j.coord] ,y = helix[k,ent.i.y],flip = flip,col = helix[k,col],lty =lty,shape=shape,lwd = lwd)
        }
      }
    }else{
      for (k in 1:helix[,.N]) {
        plotArc(i = helix[k,i.coord],j =helix[k,j.coord] ,y = y,flip = flip,col = helix[k,col],lty =lty,shape=shape,lwd = lwd)
      }
    }
  }else{
    if("ent.i.y" %in% colnames(helix)){
      for (k in 1:helix[,.N]) {
        if(helix[k,ent.i.y]!=helix[k,ent.j.y]){
          plotArcDouble(i = helix[k,i.coord],j =helix[k,j.coord] ,y1 = helix[k,ent.i.y],y2 = helix[k,ent.j.y],flip = flip,col = col,lty =lty,shape=shape,lwd = lwd)
        }else{
          plotArc(i = helix[k,i.coord],j =helix[k,j.coord] ,y = helix[k,ent.i.y],flip = flip,col = col,lty =lty,shape=shape,lwd = lwd)
        }
      }
    }else{
      for (k in 1:helix[,.N]) {
        plotArc(i = helix[k,i.coord],j =helix[k,j.coord] ,y = y,flip = flip,col = col,lty =lty,shape=shape,lwd = lwd)
      }
    }
  }
}


###############################################################################
#prep.uniq.single.reverse
###############################################################################
prep.uniq.single.reverse <- function(uniq,dist.part,full.len = 0){
  #set index N for location in plot
  index <- uniq[,.I]
  uniq[,N:=index]
  len <- uniq[,as.numeric(width)]
  uniq[,width:=len]
  #distance between entities
  dist.between <- as.integer(dist.part*mean(uniq[,width]))
  #define start and end pos
  if(full.len == 0){
    full.len <- sum(uniq[,width]) + (uniq[,.N]-1)*dist.between
  }
  for (i in 1:uniq[,.N]) {
    uniq[i,start := (uniq[,.N]-i)*dist.between+
           sum(uniq$width[i:uniq[,.N]])]
  }
  uniq[,end := start-width]
  attr(uniq,"full.length") <- full.len
  attr(uniq,"dist.between") <- dist.between
  return(uniq)
}

###############################################################################
#update.helix.coord.single.reverse
###############################################################################
update.helix.coord.single.reverse <- function(helix,uniq){
  for(k in 1:uniq[,.N]) {
    helix[ent.i %like% uniq[k,group],i.coord := uniq[k,start]-i]
    helix[ent.j %like% uniq[k,group],j.coord := uniq[k,start]-j]
  }
  return(helix)
}


###############################################################################
#update.helix.coord.single.reverse
###############################################################################
plot.arcs.contact <- function(helix,col = "black",flip = FALSE,y = 0,lty=1,shape="triangle"){

  for (k in 1:helix[,.N]){
    i <- as.numeric(helix[k,i.coord])
    j <- as.numeric(helix[k,j.coord])
    y1 <- as.numeric(helix[k,ent.i.y])
    y2 <- as.numeric(helix[k,ent.j.y])

    if(shape == "triangle"){
      if(flip){
        if("col" %in% names(helix)){
          lines(c(i, i, j), c(y1, y2, y2), col = helix$col[k])
        }else{
          lines(c(i, i, j), c(y1, y2, y2), col = col)
        }
      }else{
        if("col" %in% names(helix)){
          lines(c(i, j, j), c(y1, y1, y2), col = helix$col[k])
        }else{
          lines(c(i, j, j), c(y1, y1, y2), col = col)
        }
      }
    }


    if(shape == "heatmap"){
      if(flip){
        if("col" %in% names(helix)){
          rect(xleft = i-0.5,xright = i+0.5,ybottom = y2-0.5,ytop = y2+0.5,border = NA,col = helix$col[k])
        }else{
          rect(xleft = i-0.5,xright = i+0.5,ybottom = y2-0.5,ytop = y2+0.5,border = NA,col = col)
        }
      }else{
        if("col" %in% names(helix)){
          rect(xleft = j-0.5,xright = j+0.5,ybottom = y1-0.5,ytop = y1+0.5,border = NA,col = helix$col[k])
        }else{
          rect(xleft = j-0.5,xright = j+0.5,ybottom = y1-0.5,ytop = y1+0.5,border = NA,col = col)
        }
      }
    }
  }
}



###############################################################################
#expandNumberedHelixMltp
###############################################################################
expandNumberedHelixMltp <- function(helix) {
  number <- rep(1:nrow(helix), helix$length)
  out <- expandHelixMltp(helix)
  out$number <- number
  return(out)
}





###############################################################################
#plotScaleMod
###############################################################################
plotScaleMod <- function(width, top, bottom, ...) {
  pt <- strheight("0")
  bars <- pretty(c(1, width), n = 10)
  bars <- bars[which(bars <= width)]
  segments(bars, bottom - 1, y1 = top + 1.5 * pt, ...)
  #text(bars, top + 1.5 * pt, bars, adj = c(-0.25, 1), ...)
}


###############################################################################
#blankPlotMod
###############################################################################
#' @export
#'
blankPlotMod <- function(width, top, bottom, pad = c(0, 0, 0, 0),
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
      if(x > 2000 | y >2000){
        pdf(pdf, x/as.integer(x/1000) * factor, y/as.integer(y/1000) * factor)
      }else{
        pdf(pdf, x * factor, y * factor)
      }
    }
  }
  if (!no.par) { par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), ...) }
  plot(c(0 - pad[2], width + pad[4]), c(bottom - pad[1], top + pad[3]), type = "n", axes = F,
       xlab = "", ylab = "", asp = 1)
  if (scale) {
    plotScaleMod(width, top, bottom, col = scale.col, lwd = scale.lwd, cex = scale.cex, xpd = TRUE)
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




###############################################################################
#max.height.double
###############################################################################
max.height.double <- function(helix, msa,comp.seq = FALSE,y = 0, msa.all.seq = FALSE,top.name = "FALSE",
                              sort = FALSE, conflict.cutoff = 0.01, stable = FALSE,dist.between = 0.1,
                              dist.part = 0.2,pad = c(0, 0, 0, 0), species = 10,arcs = TRUE, palette = NA,
                              base.colour = FALSE,grid = TRUE){

  if(missing(msa)){

    options(warn = -1)
    if (helix[length > 1, .N] != 0) {
      stop("Helix is not expanded")
    }
    if (is.null(attr(helix, "length"))) {
      stop("No attribute for helix")
    }
    uniq <- attr(helix, "length")
    uniq <- update.uniq.by.group.single(uniq, top.name, sort)
    y.dist.between <- as.integer(sum(uniq[, as.numeric(width)]) *
                                   dist.between)
    uniq[1, `:=`(y, -y.dist.between)]
    uniq[!1, `:=`(y, rep(y.dist.between, uniq[2:uniq[,
                                                     .N], .N]))]
    uniq.part1 <- prep.uniq.single(uniq = uniq[1], dist.part = dist.part)
    dist.between.neg <- attr(uniq.part1, "dist.between")
    uniq.part2 <- prep.uniq.single(uniq = uniq[!1], dist.part = dist.part)
    full.len.neg <- uniq.part1[1, end]
    full.len.pos <- uniq.part2[uniq.part2[, .N], end]
    center <- as.integer(max(c(full.len.neg, full.len.pos))/2)
    if (full.len.neg > full.len.pos) {
      uniq.part2 <- prep.uniq.single(uniq = uniq.part2, dist.part = dist.part,
                                     full.len = center + (full.len.pos/2))
      dist.between <- attr(uniq.part1, "dist.between")
    }
    if (full.len.neg < full.len.pos) {
      uniq.part1 <- prep.uniq.single(uniq = uniq.part1, dist.part = dist.part,
                                     full.len = center + (full.len.neg/2))
      dist.between <- attr(uniq.part2, "dist.between")
    }
    uniq <- rbind(uniq.part1, uniq.part2)
    uniq[, `:=`(N, 1:uniq[, .N])]
    helix.pos <- helix[ent.i != uniq[1, group] & ent.j != uniq[1,group]]
    helix.pos <- copy(update.helix.coord.double(helix.pos, uniq,  2, uniq[, .N], mode = 2))
    helix.neg <- helix[ent.i == uniq[1, group] & ent.j == uniq[1,  group]]
    helix.neg <- copy(update.helix.coord.double(helix.neg, uniq,  1, 1, mode = 1))
    helix.trans <- helix[ent.i != ent.j][ent.i == uniq[1, group] | ent.j == uniq[1, group]]
    helix.trans <- copy(update.helix.coord.double(helix.trans,  uniq, 1, uniq[, .N], mode = 3))
    if (length(y) != 0) {
      if (length(y) == 1) {
        y.update <- rep(y, nrow(uniq))
        uniq[, `:=`(y.arcs, y.update)]
      }else {
        if (nrow(uniq) == length(y)) {
          y.update <- y
          uniq[, `:=`(y.arcs, y.update)]
        }else {
          stop("vector of y is incompatible with number of entities")
        }
      }
      for (i in 1:uniq[, .N]) {
        if (uniq[i, y] < 0) {
          uniq[i, `:=`(y, y - y.arcs)]
        }else {
          uniq[i, `:=`(y, y + y.arcs)]
        }
      }
      for (k in 1:uniq[, .N]) {
        if (k == 1) {
          tmp <- uniq[k, y]
          helix.neg[, `:=`(ent.i.y, tmp)]
          helix.neg[, `:=`(ent.j.y, tmp)]
        }else {
          tmp <- uniq[k, y]
          helix.pos[ent.i %like% uniq[k, group], `:=`(ent.i.y, tmp)]
          helix.pos[ent.j %like% uniq[k, group], `:=`(ent.j.y, tmp)]
        }
      }
    }
    if (stable) {
      if (nrow(helix.pos) == 0) {
        value.to.add.pos = y.dist.between
      }else {
        value.to.add.pos = max(helix.pos[, c(ent.i.y, ent.j.y)])
      }
      if (nrow(helix.neg) == 0) {
        value.to.add.neg = -y.dist.between
      }else {
        value.to.add.neg = min(helix.neg[, c(ent.i.y, ent.j.y)])
      }
      pos.max.height <- full.len.pos/2 + value.to.add.pos
      neg.max.height <- -full.len.neg/2 + value.to.add.neg
    }else {
      if (nrow(helix.pos) == 0) {
        value.to.add.pos = full.len.pos/2 + y.dist.between
      }else {
        value.to.add.pos = max(helix.pos[, abs(j.coord - i.coord)])/2 + max(helix.pos[, c(ent.i.y, ent.j.y)])
      }
      if (nrow(helix.neg) == 0) {
        value.to.add.neg = -full.len.neg/2 - y.dist.between
      }else {
        value.to.add.neg = -max(helix.neg[, abs(j.coord - i.coord)])/2 + min(helix.neg[, c(ent.i.y, ent.j.y)])
      }
      pos.max.height <- value.to.add.pos
      neg.max.height <- value.to.add.neg
    }

    return(c(pos.max.height,neg.max.height))


  }else{

    if (!is.helix.mltp(helix)) {
      stop("invalid helix input")
    }
    if (is(msa, "XStringSet")) {
      msa <- as.character(msa)
    }
    options(warn = -1)
    msa.data.table <- bstring.to.data.table(msa, comp.seq = comp.seq)
    uniq.fake <- copy(attr(helix, "length"))
    if (length(unique(uniq.fake[, group])) != length(unique(msa.data.table[,
                                                                           group]))) {
      stop("Not equal amount of entities between fasta and msa")
    }
    if (!unique(sort(unique(msa.data.table[, group])) == sort(unique(msa.data.table[,
                                                                                    group])))) {
      stop("Entities names are not same between fasta and msa")
    }
    uniq.fake <- uniq.fake[, `:=`(("width"), lapply(.SD,
                                                    as.numeric)), .SDcols = "width"]
    index <- 1:uniq.fake[, .N]
    uniq.fake[, `:=`(N, index)]
    uniq.fake <- update.uniq.by.group.single(uniq.fake, top.name,
                                             sort)
    msa.group.seq.full <- list()
    seq.names.msa <- list()
    for (i in 1:nrow(uniq.fake)) {
      msa.group.seq.full[[uniq.fake[i, group]]] <- msa.data.table[group %like%
                                                                    uniq.fake[i, group]][, !"group"][, !"rest"]
      if (comp.seq) {
        seq.names.msa[[uniq.fake[i, group]]] <- msa.group.seq.full[[uniq.fake[i,
                                                                              group]]][, seq.name]
      }
      if (!comp.seq) {
        seq.names.msa[[uniq.fake[i, group]]] <- msa.group.seq.full[[uniq.fake[i,
                                                                              group]]][, name]
      }
    }
    for (i in 1:nrow(uniq.fake)) {
      if (length(unique(duplicated(seq.names.msa[[uniq.fake[i,
                                                            group]]]))) > 1) {
        message("\nDuplicate sequence name present.\nJust first sequence will be used")
      }
      seq.names.msa[[uniq.fake[i, group]]] <- unique(seq.names.msa[[uniq.fake[i,group]]])
    }
    if (comp.seq) {
      common.seq.names <- Reduce(intersect, seq.names.msa)
      if (length(common.seq.names) == 0) {
        warning("No similar seq names present in entities\nFirst will be choosen")
      }
    }
    msa.rows.height <- c()
    for (i in 1:nrow(uniq.fake)) {
      msa.rows.height <- c(msa.rows.height, nrow(msa.group.seq.full[[uniq.fake[i,group]]]))
    }
    height.msa.min <- min(msa.rows.height)
    height.msa.max <- max(msa.rows.height)
    if (comp.seq) {
      if (height.msa.min != length(common.seq.names)) {
        message("Not enough of common sequence names between all entities\nNot common will be selected")
      }
      if (length(common.seq.names) > 0) {
        msa.selected.list <- list()
        for (i in 1:nrow(uniq.fake)) {
          tmp <- msa.group.seq.full[[uniq.fake[i, group]]][seq.name %in%
                                                             common.seq.names][!duplicated(seq.name)]
          if (nrow(tmp) != height.msa.min) {
            message("\nNot enough common sequences between entities,several will be added")
            tmp <- rbind(tmp, msa.group.seq.full[[uniq.fake[i,
                                                            group]]][!msa.selected.list[[uniq.fake[i,
                                                                                                   group]]], on = c("seq.name", "sequence")][1:(height.msa.min -
                                                                                                                                                  nrow(tmp))])
          }
          msa.selected.list[[uniq.fake[i, group]]] <- tmp[order(seq.name)]
        }
      }
      if (length(common.seq.names) == 0) {
        msa.selected.list <- list()
        for (i in 1:nrow(uniq.fake)) {
          tmp <- msa.group.seq.full[[uniq.fake[i, group]]][order(seq.name)][1:height.msa.min]
          msa.selected.list[[uniq.fake[i, group]]] <- tmp
        }
      }
    }
    else {
      msa.selected.list <- list()
      for (i in 1:nrow(uniq.fake)) {
        tmp <- msa.group.seq.full[[uniq.fake[i, group]]][order(name)][1:height.msa.min]
        msa.selected.list[[uniq.fake[i, group]]] <- tmp
      }
    }
    seq <- c()
    for (k in 1:uniq.fake[, .N]) {
      tmp <- msa.selected.list[[uniq.fake[k, group]]]$sequence
      seq <- paste(seq, tmp, sep = "")
    }
    msa.prep.equal <- BStringSet(x = seq)
    helix.trans.def <- copy(helix)
    full.len.fake <- sum(uniq.fake[, width])
    for (i in 1:uniq.fake[, .N]) {
      uniq.fake[i, `:=`(start, full.len.fake - sum(uniq.fake$width[i:uniq.fake[,
                                                                               .N]]))]
    }
    uniq.fake[, `:=`(end, width + start)]
    helix.fake <- update.helix.coord.single(helix.trans.def,
                                            uniq.fake)
    helix.fake[, `:=`(i, i.coord)]
    helix.fake[, `:=`(j, j.coord)]
    helix.fake <- helix.fake[, `:=`(c("i", "j"),
                                    lapply(.SD, as.integer)), .SDcols = c("i", "j")]
    helix.fake[, `:=`(conflict, isConflictingHelix(helix.fake[,
                                                              1:4]))]
    helix.fake.restored <- copy(helix.fake)
    helix.fake.restored[, `:=`(i, helix[, i])][, `:=`(j,
                                                      helix[, j])]
    helix.fake.restored <- helix.fake.restored[, !c("i.coord",
                                                    "j.coord")]
    helix.conflict <- helix.fake.restored[conflict >= conflict.cutoff]
    helix.notconflict <- helix.fake.restored[conflict < conflict.cutoff]
    if (base.colour) {
      cols <- getBaseColours(msa.prep.equal)
    }
    else {
      cols <- getCovarianceColoursMltp(msa.prep.equal, helix.fake[,
                                                                  1:4])
    }
    height.msa <- ifelse(grid, 1.25, 1)
    if (msa.all.seq) {
      msa.size <- height.msa * height.msa.max
    }
    else {
      msa.size <- height.msa * height.msa.min
    }
    uniq <- attr(helix, "length")
    uniq <- update.uniq.by.group.single(uniq, top.name, sort)
    uniq <- prep.uniq.single(uniq, dist.part)
    if (base.colour) {
      cols <- getBaseColours(msa.prep.equal)
    }
    else {
      cols <- getCovarianceColoursMltp(msa.prep.equal, helix.fake[,
                                                                  1:4])
    }
    if (all(is.na(palette))) {
      if (base.colour) {
        palette <- c("#E41A1C", "#4DAF4A", "#FF7F00",
                     "#377EB8", "#BDBDBD", "#984EA3")
      }
      else {
        palette <- c("#00A651", "#0072BC", "#00B9F2",
                     "#F15A22", "#231F20", "#AAAAAA",
                     "#DA6FAB")
      }
    }
    else {
      if (length(palette) < ifelse(base.colour, 5, 7)) {
        palette <- rep(palette, length.out = ifelse(base.colour,
                                                    5, 7))
      }
    }
    used <- sort(unique(c(cols)))
    dim <- dim(cols)
    cols <- palette[cols]
    dim(cols) <- dim
    uniq.fake[, `:=`(start, start + 1)]
    msa.cols.storage <- list()
    msa.name.storage <- list()
    msa.seq.storage <- list()
    for (i in 1:uniq.fake[, .N]) {
      if (i == 1) {
        if (is.vector(cols[, uniq.fake[i, start]:uniq.fake[i,end]])) {
          tmp <- matrix(cols[, uniq.fake[i, start]:uniq.fake[i, end]], nrow = 1)
        }
        else {
          tmp <- cols[, uniq.fake[i, start]:uniq.fake[i,
                                                      end]]
        }
        if (nrow(msa.group.seq.full[[uniq.fake[i, group]]]) != nrow(tmp)) {
          for (j in 1:(nrow(msa.group.seq.full[[uniq.fake[i, group]]]) - nrow(tmp))) {
            tmp <- rbind(tmp, rep("#231F20", ncol(tmp)))
          }
        }
        msa.cols.storage[[i]] <- tmp
        rm(tmp)
        if (comp.seq) {
          tmp <- rbind(msa.selected.list[[uniq.fake[i,group]]], msa.group.seq.full[[uniq.fake[i, group]]]
                       [!msa.selected.list[[uniq.fake[i, group]]],
                         on = c("seq.name", "sequence")])
        }
        else {
          tmp <- rbind(msa.selected.list[[uniq.fake[i, group]]],
                       msa.group.seq.full[[uniq.fake[i,  group]]]
                       [!msa.selected.list[[uniq.fake[i, group]]],
                         on = c("name", "sequence")])
        }
        tmp1 <- strsplit(tmp[, sequence], "")
        matrix <- matrix(ncol = length(tmp1[[1]]), nrow = length(tmp1))
        for (j in 1:length(tmp1)) {
          matrix[j, ] <- tmp1[[j]]
        }
        msa.seq.storage[[i]] <- matrix
      }
      else {
        if (is.vector(cols[, uniq.fake[i, start]:uniq.fake[i,end]])) {
          tmp <- matrix(cols[, uniq.fake[i, start]:uniq.fake[i,end]], nrow = 1)
        }
        else {
          tmp <- cols[, uniq.fake[i, start]:uniq.fake[i,end]]
        }
        if (nrow(msa.group.seq.full[[uniq.fake[i, group]]]) != nrow(tmp)) {
          for (j in 1:(nrow(msa.group.seq.full[[uniq.fake[i,group]]]) - nrow(tmp))) {
            tmp <- rbind(tmp, rep("#231F20", ncol(tmp)))
          }
        }
        for (j in 1:nrow(tmp)) {
          tmp[j, ] <- rev(tmp[j, ])
        }
        msa.cols.storage[[i]] <- tmp
        if (comp.seq) {
          tmp <- rbind(msa.selected.list[[uniq.fake[i,group]]],
                       msa.group.seq.full[[uniq.fake[i,group]]]
                       [!msa.selected.list[[uniq.fake[i, group]]],
                         on = c("seq.name", "sequence")])
        }
        else {
          tmp <- rbind(msa.selected.list[[uniq.fake[i,group]]],
                       msa.group.seq.full[[uniq.fake[i,group]]]
                       [!msa.selected.list[[uniq.fake[i, group]]],
                         on = c("name", "sequence")])
        }
        tmp1 <- strsplit(tmp[, sequence], "")
        matrix <- matrix(ncol = length(tmp1[[1]]), nrow = length(tmp1))
        for (j in 1:length(tmp1)) {
          matrix[j, ] <- rev(tmp1[[j]])
        }
        msa.seq.storage[[i]] <- matrix
      }
      if (comp.seq) {
        msa.name.storage[[i]] <- tmp[, seq.name]
      }
      else {
        msa.name.storage[[i]] <- tmp[, name]
      }
    }
    if (!msa.all.seq) {
      for (i in 1:length(msa.cols.storage)) {
        msa.seq.storage[[i]] <- msa.seq.storage[[i]][1:height.msa.min,
                                                     ]
        msa.cols.storage[[i]] <- msa.cols.storage[[i]][1:height.msa.min,
                                                       ]
        msa.name.storage[[i]] <- msa.name.storage[[i]][1:height.msa.min]
      }
    }
    y.msa <- c()
    for (i in 1:length(msa.cols.storage)) {
      if (is.vector(msa.cols.storage[[i]])) {
        y.msa <- c(y.msa, 1 * height.msa)
      }
      else {
        y.msa <- c(y.msa, nrow(msa.cols.storage[[i]]) * height.msa)
      }
    }
    uniq <- copy(attr(helix, "length"))
    uniq <- update.uniq.by.group.single(uniq, top.name, sort)
    y.dist.between <- as.integer(sum(uniq[, as.numeric(width)]) *
                                   dist.between)
    uniq[1, `:=`(y, -y.dist.between)]
    uniq[!1, `:=`(y, rep(y.dist.between, uniq[2:uniq[,
                                                     .N], .N]))]
    uniq.part1 <- prep.uniq.single(uniq = uniq[1], dist.part = dist.part)
    dist.between.neg <- attr(uniq.part1, "dist.between")
    uniq.part2 <- prep.uniq.single(uniq = uniq[!1], dist.part = dist.part)
    full.len.neg <- uniq.part1[1, end]
    full.len.pos <- uniq.part2[uniq.part2[, .N], end]
    center <- as.integer(max(c(full.len.neg, full.len.pos))/2)
    if (full.len.neg > full.len.pos) {
      uniq.part2 <- prep.uniq.single(uniq = uniq.part2, dist.part = dist.part,
                                     full.len = center + (full.len.pos/2))
      dist.between.uniq <- attr(uniq.part1, "dist.between")
    }
    if (full.len.neg < full.len.pos) {
      uniq.part1 <- prep.uniq.single(uniq = uniq.part1, dist.part = dist.part,
                                     full.len = center + (full.len.neg/2))
      dist.between.uniq <- attr(uniq.part2, "dist.between")
    }
    uniq <- rbind(uniq.part1, uniq.part2)
    uniq[, `:=`(N, 1:uniq[, .N])]
    helix.pos <- helix[ent.i != uniq[1, group] & ent.j != uniq[1,
                                                               group]]
    helix.pos <- update.helix.coord.double(helix.pos, uniq, 2,
                                           uniq[, .N], mode = 2)
    helix.neg <- helix[ent.i == uniq[1, group] & ent.j == uniq[1,
                                                               group]]
    helix.neg <- update.helix.coord.double(helix.neg, uniq, 1,
                                           1, mode = 1)
    helix.trans <- helix[ent.i != ent.j][ent.i == uniq[1, group] |
                                           ent.j == uniq[1, group]]
    helix.trans <- update.helix.coord.double(helix.trans, uniq,
                                             1, uniq[, .N], mode = 3)
    if (stable) {
      pos.max.height <- full.len.pos/2 + y.dist.between
      neg.max.height <- -full.len.neg/2 - y.dist.between
    }else {
      if (nrow(helix.pos) == 0) {
        value.to.add.pos = full.len.pos/2 + y.dist.between
      }else {
        value.to.add.pos = max(helix.pos[, abs(j.coord - i.coord)])/2 + y.dist.between
      }
      if (nrow(helix.neg) == 0) {
        value.to.add.neg = -full.len.neg/2 - y.dist.between
      }else {
        value.to.add.neg = -max(helix.neg[, abs(j.coord - i.coord)])/2 - y.dist.between
      }
      pos.max.height <- value.to.add.pos
      neg.max.height <- value.to.add.neg
    }
    height.msa <- ifelse(grid, 1.25, 1)
    if (msa.all.seq) {
      neg.msa <- nrow(msa.seq.storage[[1]])
      tmp <- c()
      for (i in 2:length(msa.seq.storage)) {
        tmp <- c(tmp, nrow(msa.seq.storage[[i]]))
      }
      pos.msa <- c(max(tmp))
      if (arcs) {
        pos.max.height <- pos.max.height + height.msa * pos.msa
        neg.max.height <- neg.max.height - height.msa * neg.msa
      }else {
        pos.max.height <- y.dist.between + height.msa * pos.msa
        neg.max.height <- -y.dist.between - height.msa *
          neg.msa
      }
    }else {
      if (arcs) {
        pos.max.height <- pos.max.height + height.msa * height.msa.min
        neg.max.height <- neg.max.height - height.msa * height.msa.min
      }
      else {
        pos.max.height <- y.dist.between + height.msa * height.msa.min
        neg.max.height <- -y.dist.between - height.msa * height.msa.min
      }
    }


    return(c(pos.max.height,neg.max.height))
  }


}



###############################################################################
#max.width.double
###############################################################################
max.width.double <- function(helix, msa,comp.seq = FALSE,y = 0, msa.all.seq = FALSE,top.name = "FALSE",
                             sort = FALSE, conflict.cutoff = 0.01, stable = FALSE,dist.between = 0.1,
                             dist.part = 0.2,pad = c(0, 0, 0, 0), species = 10,arcs = TRUE, palette = NA){

  if(missing(msa)){

    options(warn = -1)
    if (helix[length > 1, .N] != 0) {
      stop("Helix is not expanded")
    }
    if (is.null(attr(helix, "length"))) {
      stop("No attribute for helix")
    }
    uniq <- attr(helix, "length")
    uniq <- update.uniq.by.group.single(uniq, top.name, sort)
    y.dist.between <- as.integer(sum(uniq[, as.numeric(width)]) *
                                   dist.between)
    uniq[1, `:=`(y, -y.dist.between)]
    uniq[!1, `:=`(y, rep(y.dist.between, uniq[2:uniq[,
                                                     .N], .N]))]
    uniq.part1 <- prep.uniq.single(uniq = uniq[1], dist.part = dist.part)
    dist.between.neg <- attr(uniq.part1, "dist.between")
    uniq.part2 <- prep.uniq.single(uniq = uniq[!1], dist.part = dist.part)
    full.len.neg <- uniq.part1[1, end]
    full.len.pos <- uniq.part2[uniq.part2[, .N], end]
    center <- as.integer(max(c(full.len.neg, full.len.pos))/2)
    if (full.len.neg > full.len.pos) {
      uniq.part2 <- prep.uniq.single(uniq = uniq.part2, dist.part = dist.part,
                                     full.len = center + (full.len.pos/2))
      dist.between <- attr(uniq.part1, "dist.between")
    }
    if (full.len.neg < full.len.pos) {
      uniq.part1 <- prep.uniq.single(uniq = uniq.part1, dist.part = dist.part,
                                     full.len = center + (full.len.neg/2))
      dist.between <- attr(uniq.part2, "dist.between")
    }

    return(max(full.len.neg,full.len.pos))


  }else{

    uniq <- copy(attr(helix, "length"))
    uniq <- update.uniq.by.group.single(uniq, top.name, sort)
    y.dist.between <- as.integer(sum(uniq[, as.numeric(width)]) *
                                   dist.between)
    uniq[1, `:=`(y, -y.dist.between)]
    uniq[!1, `:=`(y, rep(y.dist.between, uniq[2:uniq[,
                                                     .N], .N]))]
    uniq.part1 <- prep.uniq.single(uniq = uniq[1], dist.part = dist.part)
    dist.between.neg <- attr(uniq.part1, "dist.between")
    uniq.part2 <- prep.uniq.single(uniq = uniq[!1], dist.part = dist.part)
    full.len.neg <- uniq.part1[1, end]
    full.len.pos <- uniq.part2[uniq.part2[, .N], end]
    center <- as.integer(max(c(full.len.neg, full.len.pos))/2)
    if (full.len.neg > full.len.pos) {
      uniq.part2 <- prep.uniq.single(uniq = uniq.part2, dist.part = dist.part,
                                     full.len = center + (full.len.pos/2))
      dist.between.uniq <- attr(uniq.part1, "dist.between")
    }
    if (full.len.neg < full.len.pos) {
      uniq.part1 <- prep.uniq.single(uniq = uniq.part1, dist.part = dist.part,
                                     full.len = center + (full.len.neg/2))
      dist.between.uniq <- attr(uniq.part2, "dist.between")
    }

    return(max(full.len.neg,full.len.pos))

  }
}
