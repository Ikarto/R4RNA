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








#' @rdname plotHelixMltpSingleLine
#' @export
#'
###############################################################################
#plotHelixMltpDoubleLine
###############################################################################
plotHelixMltpDoubleLine <- function(helix,top.name = "FALSE",  sort = FALSE,
                                    dist.between = .1,  dist.part = .2,
                                    stable = FALSE,  line = TRUE,
                                    arrow = TRUE,  col.line = "black",
                                    col.arrow = "black",  col = "black",
                                    shape = "circle",  scale = FALSE,
                                    debug = FALSE,add = FALSE,y=0,arc.lty=1,append = FALSE,x = 0,...){

  #restore warnings
  options(warn = -1)
  if (helix[length > 1, .N] != 0) {
    stop("Helix is not expanded")
  }
  if (is.null(attr(helix, "length"))) {
    stop("No attribute for helix")
  }
  uniq <- attr(helix, "length")
  uniq <- update.uniq.by.group.single(uniq, top.name, sort)
  y.dist.between <- as.integer(sum(uniq[, as.numeric(width)]) * dist.between)
  uniq[1, `:=`(y, -y.dist.between)]
  uniq[!1, `:=`(y, rep(y.dist.between, uniq[2:uniq[,.N], .N]))]
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
  helix.pos <- copy(update.helix.coord.double(helix.pos, uniq, 2, uniq[, .N], mode = 2))
  helix.neg <- helix[ent.i == uniq[1, group] & ent.j == uniq[1, group]]
  helix.neg <- copy(update.helix.coord.double(helix.neg, uniq,1, 1, mode = 1))
  helix.trans <- helix[ent.i != ent.j][ent.i == uniq[1, group] | ent.j == uniq[1, group]]
  helix.trans <- copy(update.helix.coord.double(helix.trans, uniq, 1, uniq[, .N], mode = 3))

  if (y != 0) {
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
      }
      else {
        tmp <- uniq[k, y]
        helix.pos[ent.i %like% uniq[k, group], `:=`(ent.i.y, tmp)]
        helix.pos[ent.j %like% uniq[k, group], `:=`(ent.j.y, tmp)]
      }
    }
  }


  if (stable) {
    if(nrow(helix.pos) == 0){
      value.to.add.pos = y.dist.between
    }else{
      value.to.add.pos = max(helix.pos[, c(ent.i.y, ent.j.y)])
    }

    if(nrow(helix.neg) == 0){
      value.to.add.neg = -y.dist.between
    }else{
      value.to.add.neg = min(helix.neg[, c(ent.i.y, ent.j.y)])
    }


    pos.max.height <- full.len.pos/2 + value.to.add.pos
    neg.max.height <- -full.len.neg/2 + value.to.add.neg
  }else {
    if(nrow(helix.pos) == 0){
      value.to.add.pos = full.len.pos/2 + y.dist.between
    }else{
      value.to.add.pos = max(helix.pos[, abs(j.coord - i.coord)])/2 + max(helix.pos[, c(ent.i.y, ent.j.y)])
    }

    if(nrow(helix.neg) == 0){
      value.to.add.neg = -full.len.neg/2 - y.dist.between
    }else{
      value.to.add.neg = -max(helix.neg[, abs(j.coord - i.coord)])/2 + min(helix.neg[, c(ent.i.y, ent.j.y)])
    }


    pos.max.height <- value.to.add.pos
    neg.max.height <- value.to.add.neg
  }



  if(x != 0){
    uniq <- uniq[,start := start + x]
    uniq <- uniq[,end := end + x]


    if (helix.pos[, .N] > 0) {
      helix.pos <- helix.pos[, i.coord := i.coord + x ]
      helix.pos <- helix.pos[, j.coord := j.coord + x ]
    }
    if (helix.neg[, .N] > 0) {
      helix.neg <- helix.neg[, i.coord := i.coord + x ]
      helix.neg <- helix.neg[, j.coord := j.coord + x ]
    }
    if (helix.trans[, .N] > 0) {
      helix.trans <- helix.trans[, i.coord := i.coord + x ]
      helix.trans <- helix.trans[, j.coord := j.coord + x ]
    }

  }


  if (!add) {
    if (append) {
      blankPlotMod(max(c(full.len.neg, full.len.pos)) + dist.between/8 + x,
                   pos.max.height, neg.max.height, scale = scale,
                   debug = debug, no.par = TRUE,...)
    }else {
      blankPlotMod(max(c(full.len.neg, full.len.pos)) + dist.between/8 + x,
                   pos.max.height, neg.max.height, scale = scale,
                   debug = debug,...)
    }
  }
  if (line) {
    plot.bases.single(uniq, col.line)
  }
  if (arrow) {
    plot.arrows.single(uniq, y.dist.between, col.arrow, mode = 2)
  }
  if (helix.pos[, .N] > 0) {
    plot.arcs.single(helix.pos, col = col, flip = FALSE,
                     lty = arc.lty, shape = shape, ...)
  }
  if (helix.neg[, .N] > 0) {
    plot.arcs.single(helix.neg, col = col, flip = TRUE, lty = arc.lty,
                     shape = shape,...)
  }
  if (helix.trans[, .N] > 0) {
    plot.lines.single(helix.trans, col = col, lty = arc.lty,...)
  }
  #restore warnings
  options(warn = 1)
}





#' Plot Helix in arc diagram in double line mode in double mode
#'
#' @description
#'
#' Plots 2 helix data tables as an arc diagrams in double line (one entity against others,
#' where one entity can be choosen by top.name option),
#' with styling possible with properly named additional columns on the data table separately.
#'
#' @param left,right helixs data.table, with the six mandatory columns. "col" column can be specified as styling column.
#' @param top.name Name of entity, which will be set as top name and will refer to first entity in plot
#' @param sort sort entities by alphabetic order
#' @param dist.part distance between entities, dist.part*average length of entities
#' @param dist.between distance between entities on y axes, dist.between*average length of entities
#' @param stable used to have same height of resulted plot, based on length of entities
#' @param line If TRUE, a horizontal line representing the sequence is plotted.
#' @param arrow If TRUE, an arrow is played on the right end of the line.
#' @param col.line Colour of line, default is "black"
#' @param col.arrow Colour of arrow, default is "black"
#' @param col colour of arcs, if there is no style column in helix data.table sstructure
#' @param shape One of "circle", "triangle", or "square", specifying the shape of the arcs.
#' @param scale If TRUE, inserts a scale on the plot.
#' @param debug If TRUE, frames the boundaries of the intended plotting space in red, used to determine if inputs produce expected output area. Also outputs to STDIN dimensions of the plot.
#' @param append If TRUE, graphical elements are added to the active plot device, else a new plot device is created for the plot.
#'
#' @export
#'
###############################################################################
#plotDoubleHelixMltpDoubleLine
###############################################################################
plotDoubleHelixMltpDoubleLine <- function(left,right,top.name = "FALSE",
                                          sort = TRUE,dist.between = .1,
                                          dist.part = .2,dist.par.plot = .2,stable = FALSE,
                                          line = TRUE,arrow = TRUE,
                                          col.line = "black",col.arrow = "black",
                                          col = "black",shape = "circle",
                                          scale = FALSE,debug = FALSE,append=TRUE,...){



  height1 <- max.height.double(helix = left, comp.seq = comp.seq,top.name = top.name,
                               sort = sort, stable = stable,dist.between = dist.between,
                               dist.part = dist.part)

  height2 <- max.height.double(helix = right, comp.seq = comp.seq,top.name = top.name,
                               sort = sort, stable = stable,dist.between = dist.between,
                               dist.part = dist.part)
  width1 <- max.width.double(helix = left,comp.seq = comp.seq,top.name = top.name,sort = sort,
                             stable = stable,dist.between = dist.between, dist.part = dist.part)

  width2 <- max.width.double(helix = right,comp.seq = comp.seq,top.name = top.name,sort = sort,
                             stable = stable,dist.between = dist.between, dist.part = dist.part)


  blankPlotMod(width = width1+width2+(max(width1,width2)*dist.par.plot),
               top = max(height1[1],height2[1]),bottom = min(height1[2],height2[2]),
               scale = scale,debug = debug,...)

  plotHelixMltpDoubleLine(left,top.name = top.name,sort = sort,
                          dist.between = dist.between,dist.part = dist.part,
                          stable = stable,  line = line,
                          arrow = arrow,  col.line = col.line,
                          col.arrow = col.arrow,  col = col,
                          shape = shape,  scale = scale,
                          debug = debug,add=TRUE,append = FALSE)

  dist.part = dist.part
  dist.between = dist.between

  plotHelixMltpDoubleLine(right,top.name = top.name,sort = sort,
                          dist.between = dist.between,dist.part = dist.part,
                          stable = stable,  line = line,
                          arrow = arrow,  col.line = col.line,
                          col.arrow = col.arrow,  col = col,
                          shape = shape,  scale = scale,
                          debug = debug,add=TRUE,append = FALSE,
                          x = width1 + (max(width1,width2)*dist.par.plot))

}


#' @rdname plotCovarianceMltpSingleLine
#' @export
#'
###############################################################################
#plotCovarianceMltpDoubleLine
###############################################################################
plotCovarianceMltpDoubleLine <- function(helix,msa,comp.seq = FALSE,
                                         y=0,msa.all.seq = FALSE,arcs = TRUE,
                                         top.name = "FALSE",sort = FALSE,
                                         conflict.cutoff = 0.01,stable = FALSE,grid = TRUE,
                                         dist.between = .1,dist.part = .2,
                                         add = FALSE,append = FALSE,pad = c(0, 0, 0, 0),
                                         shape = "circle",debug = FALSE,
                                         col = "black",conflict.col = NA,
                                         conflict.lty = 2,conflict.lwd = 1,base.colour = FALSE,
                                         palette =NA,text = TRUE, grid.col = "white",
                                         grid.lwd = 0,text.cex = 0.2,text.col = "gray",
                                         text.font = 2,text.family = "sans",
                                         species.cex = 0.2,species.col = "black",
                                         species.font = 2,species.family = "mono",
                                         species = 10,legend = TRUE,col.line = "black",
                                         col.arrow = "black",x = 0,...){

  if(!is.helix.mltp(helix)){
    stop("invalid helix input")
  }
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }
  #restore warnings
  options(warn = -1)

  #####Extract sequences and make msa as 1 BStringSet
  #####For next using on default functions
  #transfor fasta to data.table format
  msa.data.table <- bstring.to.data.table(msa,comp.seq = comp.seq)
  #call uniq from helix file
  uniq.fake <- copy(attr(helix,"length"))

  #check if N of unique entities are same
  if(length(unique(uniq.fake[,group])) != length(unique(msa.data.table[,group]))){
    stop("Not equal amount of entities between fasta and msa")
  }
  if(!unique(sort(unique(msa.data.table[,group])) == sort(unique(msa.data.table[,group])))){
    stop("Entities names are not same between fasta and msa")
  }

  uniq.fake <- uniq.fake[,("width"):= lapply(.SD, as.numeric), .SDcols = "width"]
  index <- 1:uniq.fake[,.N]
  uniq.fake[,N:=index]
  #update uniq if it needed by top.name
  uniq.fake <- update.uniq.by.group.single(uniq.fake,top.name,sort)

  #contain full msa sequences
  msa.group.seq.full <- list()
  #contain names of sequences
  seq.names.msa <- list()
  for(i in 1:nrow(uniq.fake)){
    msa.group.seq.full[[uniq.fake[i,group]]] <- msa.data.table[group %like% uniq.fake[i,group]][,!"group"][,!"rest"]
    if(comp.seq){
      seq.names.msa[[uniq.fake[i,group]]] <- msa.group.seq.full[[uniq.fake[i,group]]][,seq.name]
    }
    if(!comp.seq){
      seq.names.msa[[uniq.fake[i,group]]] <- msa.group.seq.full[[uniq.fake[i,group]]][,name]
    }

  }
  #define duplicates and take massage to user, after delete duplicate names
  for(i in 1:nrow(uniq.fake)){
    if(length(unique(duplicated(seq.names.msa[[uniq.fake[i,group]]])))>1){
      message("\nDuplicate sequence name present.\nJust first sequence will be used")
    }
    seq.names.msa[[uniq.fake[i,group]]] <- unique(seq.names.msa[[uniq.fake[i,group]]])
  }

  if(comp.seq){
    #find common names between all entities
    common.seq.names <- Reduce(intersect, seq.names.msa)
    if(length(common.seq.names)==0){
      warning("No similar seq names present in entities\nFirst will be choosen")
    }
  }
  #define minimum of msa seq
  msa.rows.height <- c()
  for(i in 1:nrow(uniq.fake)){
    msa.rows.height <-c(msa.rows.height,nrow(msa.group.seq.full[[uniq.fake[i,group]]]))
  }
  height.msa.min <- min(msa.rows.height)
  height.msa.max <- max(msa.rows.height)

  if(comp.seq){
    if(height.msa.min != length(common.seq.names)){
      message("Not enough of common sequence names between all entities\nNot common will be selected")
    }
    if(length(common.seq.names)>0){
      #joining msa's to create same levels of comparing for all of them
      msa.selected.list <- list()
      for(i in 1:nrow(uniq.fake)){
        tmp <- msa.group.seq.full[[uniq.fake[i,group]]][seq.name %in% common.seq.names][!duplicated(seq.name)]
        if(nrow(tmp)!= height.msa.min){
          message("\nNot enough common sequences between entities,several will be added")
          tmp <- rbind(tmp,msa.group.seq.full[[uniq.fake[i,group]]][!msa.selected.list[[uniq.fake[i,group]]],on = c("seq.name","sequence")][1:(height.msa.min-nrow(tmp))])
        }
        #sort by seq.name
        msa.selected.list[[uniq.fake[i,group]]] <- tmp[order(seq.name)]
      }
    }
    if(length(common.seq.names)==0){
      msa.selected.list <- list()
      for(i in 1:nrow(uniq.fake)){
        tmp <- msa.group.seq.full[[uniq.fake[i,group]]][order(seq.name)][1:height.msa.min]
        #sort by seq.name
        msa.selected.list[[uniq.fake[i,group]]] <- tmp
      }
    }
  }else{
    msa.selected.list <- list()
    for(i in 1:nrow(uniq.fake)){
      tmp <- msa.group.seq.full[[uniq.fake[i,group]]][order(name)][1:height.msa.min]
      #sort by seq.name
      msa.selected.list[[uniq.fake[i,group]]] <- tmp
    }
  }

  #JOIN selected sequences
  seq <- c()
  for (k in 1:uniq.fake[,.N]) {
    tmp <- msa.selected.list[[uniq.fake[k,group]]]$sequence
    seq <- paste(seq,tmp,sep = "")
  }
  msa.prep.equal <- BStringSet(x = seq)

  #####Transform helix
  #####For next using on default functions
  helix.trans.def <- copy(helix)
  #update uniq with fake start and end
  full.len.fake <- sum(uniq.fake[,width])
  for (i in 1:uniq.fake[,.N]) {
    uniq.fake[i,start := full.len.fake-sum(uniq.fake$width[i:uniq.fake[,.N]])]
  }
  uniq.fake[,end := width+start]
  helix.fake <- update.helix.coord.single(helix.trans.def,uniq.fake)
  helix.fake[,i:=i.coord]
  helix.fake[,j:=j.coord]
  helix.fake <- helix.fake[,c("i","j"):= lapply(.SD, as.integer), .SDcols = c("i","j")]
  #define conflicting
  helix.fake[,conflict:=isConflictingHelix(helix.fake[,1:4])]
  #create 2 helix with conflicting and non-conflicting
  helix.fake.restored <- copy(helix.fake)
  helix.fake.restored[,i:=helix[,i]][,j:=helix[,j]]
  helix.fake.restored <-helix.fake.restored[,!c("i.coord","j.coord")]
  helix.conflict <- helix.fake.restored[conflict >= conflict.cutoff]
  helix.notconflict <- helix.fake.restored[conflict < conflict.cutoff]

  ##################msa colouring and storaging##################
  if(base.colour){
    cols <- getBaseColoursMod(msa.prep.equal)
  }else{
    cols <- getCovarianceColoursMltp(msa.prep.equal,helix.fake[,1:4])
  }

  #get height of msa
  height.msa <- ifelse(grid,1.25,1)
  #chose height by amount of plotted sequences
  if(msa.all.seq){
    msa.size <- height.msa*height.msa.max
  }else{
    msa.size <- height.msa*height.msa.min
  }

  #get start and end positions
  #options(warn = -1)
  #read len data
  uniq <- attr(helix,"length")
  #sort and define top name
  uniq <- update.uniq.by.group.single(uniq,top.name,sort)
  #calculate values for uniq
  uniq <- prep.uniq.single(uniq,dist.part)


  ####plot msa
  if(base.colour){
    cols <- getBaseColoursMod(msa.prep.equal)
  }else{
    cols <- getCovarianceColoursMltp(msa.prep.equal,helix.fake[,1:4])
  }

  #pallete colours
  if (all(is.na(palette))) {
    if (base.colour) {
      palette <- c("#E41A1C", "#4DAF4A", "#FF7F00", "#377EB8", "#BDBDBD","#984EA3","#00008b")
    } else {
      palette <- c("#00A651", "#0072BC", "#00B9F2", "#F15A22", "#231F20",
                   "#AAAAAA", "#DA6FAB")
    }
  } else {
    if (length(palette) < ifelse(base.colour, 5, 7)) {
      palette <- rep(palette, length.out = ifelse(base.colour, 5, 7))
    }
  }

  #pemade cols by pallete
  used <- sort(unique(c(cols)))
  dim <- dim(cols)
  cols <- palette[cols]
  dim(cols) <- dim

  #extart separate group colors segments and create separate
  uniq.fake[,start:= start+1]
  #list with msa's
  msa.cols.storage <- list()
  msa.name.storage <- list()
  msa.seq.storage <- list()

  #msa.selected.list msa.group.seq.full
  for (i in 1:uniq.fake[,.N]) {
    if(i==1){
      ############
      #colours
      ############
      if (is.vector(cols[,uniq.fake[i,start]:uniq.fake[i,end]])) {
        tmp <- matrix(cols[,uniq.fake[i,start]:uniq.fake[i,end]],nrow = 1)
      }else{
        tmp <- cols[,uniq.fake[i,start]:uniq.fake[i,end]]
      }

      if(nrow(msa.group.seq.full[[uniq.fake[i,group]]])!=nrow(tmp)){
        for (j in 1:(nrow(msa.group.seq.full[[uniq.fake[i,group]]])-nrow(tmp))) {
          tmp <- rbind(tmp,rep("#231F20",ncol(tmp)))
        }
      }
      msa.cols.storage[[i]] <- tmp
      rm(tmp)

      if(comp.seq){
        tmp <- rbind(msa.selected.list[[uniq.fake[i,group]]],msa.group.seq.full[[uniq.fake[i,group]]][!msa.selected.list[[uniq.fake[i,group]]],on = c("seq.name","sequence")])
      }else{
        tmp <- rbind(msa.selected.list[[uniq.fake[i,group]]],msa.group.seq.full[[uniq.fake[i,group]]][!msa.selected.list[[uniq.fake[i,group]]],on = c("name","sequence")])
      }
      ############
      #sequences
      ############
      tmp1 <- strsplit(tmp[,sequence],"")
      matrix <- matrix(ncol = length(tmp1[[1]]),nrow = length(tmp1))
      for (j in 1:length(tmp1)) {
        matrix[j,] <- tmp1[[j]]
      }
      msa.seq.storage[[i]] <- matrix
    }else{
      ############
      #colours
      ############
      if (is.vector(cols[,uniq.fake[i,start]:uniq.fake[i,end]])) {
        tmp <- matrix(cols[,uniq.fake[i,start]:uniq.fake[i,end]],nrow = 1)
      }else{
        tmp <- cols[,uniq.fake[i,start]:uniq.fake[i,end]]
      }

      if(nrow(msa.group.seq.full[[uniq.fake[i,group]]])!=nrow(tmp)){
        for (j in 1:(nrow(msa.group.seq.full[[uniq.fake[i,group]]])-nrow(tmp))) {
          tmp <- rbind(tmp,rep("#231F20",ncol(tmp)))
        }
      }
      for (j in 1:nrow(tmp)) {
        tmp[j,] <- rev(tmp[j,])
      }
      msa.cols.storage[[i]] <- tmp
      ############
      #sequences
      ############
      if(comp.seq){
        tmp <- rbind(msa.selected.list[[uniq.fake[i,group]]],msa.group.seq.full[[uniq.fake[i,group]]][!msa.selected.list[[uniq.fake[i,group]]],on = c("seq.name","sequence")])
      }else{
        tmp <- rbind(msa.selected.list[[uniq.fake[i,group]]],msa.group.seq.full[[uniq.fake[i,group]]][!msa.selected.list[[uniq.fake[i,group]]],on = c("name","sequence")])
      }

      tmp1 <- strsplit(tmp[,sequence],"")
      matrix <- matrix(ncol = length(tmp1[[1]]),nrow = length(tmp1))
      for (j in 1:length(tmp1)) {
        matrix[j,] <- rev(tmp1[[j]])
      }
      msa.seq.storage[[i]] <- matrix

    }
    ############
    #names
    ############
    if(comp.seq){
      msa.name.storage[[i]] <- tmp[,seq.name]
    }else{
      msa.name.storage[[i]] <- tmp[,name]
    }

  }

  #print(msa.name.storage)
  #cut msa if it cutted version
  if(!msa.all.seq){
    for (i in 1:length(msa.cols.storage)) {
      msa.seq.storage[[i]] <- msa.seq.storage[[i]][1:height.msa.min,]
      msa.cols.storage[[i]] <- msa.cols.storage[[i]][1:height.msa.min,]
      msa.name.storage[[i]] <- msa.name.storage[[i]][1:height.msa.min]
    }
  }

  #create y alocation vector
  y.msa <- c()
  for (i in 1:length(msa.cols.storage)) {
    if (is.vector(msa.cols.storage[[i]])) {
      y.msa <- c(y.msa,1*height.msa)
    }else{
      y.msa <- c(y.msa,nrow(msa.cols.storage[[i]])*height.msa)
    }
  }

  #############################PLOTTING#######################################
  #############define max and neg height, start points and y coordinates
  uniq <- copy(attr(helix,"length"))
  #sort and define top name
  uniq <- update.uniq.by.group.single(uniq,top.name,sort)
  #give y coordinate for uniq
  y.dist.between <- as.integer(sum(uniq[,as.numeric(width)])*dist.between)
  uniq[1,y:=-y.dist.between]
  uniq[!1,y:=rep(y.dist.between,uniq[2:uniq[,.N],.N])]
  #define full len for positive part and negative
  uniq.part1 <- prep.uniq.single(uniq = uniq[1],dist.part = dist.part)
  dist.between.neg <- attr(uniq.part1,"dist.between")
  uniq.part2 <- prep.uniq.single(uniq = uniq[!1],dist.part = dist.part)
  full.len.neg <- uniq.part1[1,end]
  full.len.pos <- uniq.part2[uniq.part2[,.N],end]
  #define center point for plot
  center <- as.integer(max(c(full.len.neg,full.len.pos))/2)
  #resize our start and end for new location
  if(full.len.neg>full.len.pos){
    uniq.part2<- prep.uniq.single(uniq = uniq.part2,dist.part = dist.part,full.len = center+(full.len.pos/2))
    dist.between.uniq <- attr(uniq.part1,"dist.between")
  }
  if(full.len.neg<full.len.pos){
    uniq.part1<- prep.uniq.single(uniq = uniq.part1,dist.part = dist.part,full.len = center+(full.len.neg/2))
    dist.between.uniq <- attr(uniq.part2,"dist.between")
  }
  uniq <- rbind(uniq.part1,uniq.part2)
  uniq[,N := 1:uniq[,.N]]
  #for positive side
  helix.pos <- helix[ent.i != uniq[1,group] & ent.j!= uniq[1,group]]
  helix.pos <- update.helix.coord.double(helix.pos,uniq,2,uniq[,.N],mode=2)
  #for negative side
  helix.neg <- helix[ent.i == uniq[1,group] & ent.j== uniq[1,group]]
  helix.neg <- update.helix.coord.double(helix.neg,uniq,1,1,mode=1)
  #for trans interactions
  helix.trans <- helix[ent.i!=ent.j][ent.i == uniq[1,group] | ent.j == uniq[1,group]]
  helix.trans <- update.helix.coord.double(helix.trans,uniq,1,uniq[,.N],mode=3)


  ##########define max positive and negative side
  if(stable){
    pos.max.height <- full.len.pos/2 + y.dist.between
    neg.max.height <- - full.len.neg/2 - y.dist.between
  }else{
    if (nrow(helix.pos) == 0) {
      value.to.add.pos = full.len.pos/2 + y.dist.between
    }else{
      value.to.add.pos = max(helix.pos[,abs(j.coord-i.coord)])/2 + y.dist.between
    }

    if (nrow(helix.neg) == 0) {
      value.to.add.neg = - full.len.neg/2 - y.dist.between
    }else{
      value.to.add.neg =  -max(helix.neg[,abs(j.coord-i.coord)])/2 - y.dist.between
    }

    pos.max.height <- value.to.add.pos
    neg.max.height <- value.to.add.neg
  }


  #update height by msa height
  height.msa <- ifelse(grid,1.25,1)
  if(msa.all.seq){
    neg.msa <- nrow(msa.seq.storage[[1]])
    tmp <- c()
    for(i in 2:length(msa.seq.storage)){
      tmp <- c(tmp,nrow(msa.seq.storage[[i]]))
    }
    pos.msa <- c(max(tmp))

    if(arcs){
      pos.max.height <- pos.max.height + height.msa*pos.msa
      neg.max.height <- neg.max.height - height.msa*neg.msa
    }else{
      pos.max.height <- y.dist.between + height.msa*pos.msa
      neg.max.height <- -y.dist.between - height.msa*neg.msa
    }

  }else{
    if(arcs){
      pos.max.height <- pos.max.height + height.msa*height.msa.min
      neg.max.height <- neg.max.height - height.msa*height.msa.min
    }else{
      pos.max.height <- y.dist.between + height.msa*height.msa.min
      neg.max.height <- -y.dist.between - height.msa*height.msa.min
    }
  }

  #update by x
  if(x != 0){
    uniq <- uniq[,start := start + x]
    uniq <- uniq[,end := end + x]

  }


  #plot blank plot
  if(!add){
    if(append){
      if(species > 0){
        pad[2] <- pad[2] + species
      }
      blankPlotMod(width = max(full.len.neg,full.len.pos)+x, top = pos.max.height,bottom = neg.max.height, pad,debug = debug,no.par = TRUE,...)
    }else{
      if(species > 0){
        pad[2] <- pad[2] + species
      }
      blankPlotMod(width = max(full.len.neg,full.len.pos)+x, top = pos.max.height,bottom = neg.max.height, pad,debug = debug,...)
    }
  }


  if(arcs){
    if (!is.na(conflict.col)){
      if(!is.na(conflict.lty)) {
        if(nrow(helix.conflict) != 0){
          helix.conflict <- copy(helix.conflict[,!"conflict"])
          helix.conflict <- as.data.table(helix.conflict)
          helix.conflict <- expandHelixMltp(helix.conflict)
          attr(helix.conflict,"length") <- attr(helix,"length")
          plotHelixMltpDoubleLine(helix.conflict[,1:6],top.name = top.name,sort =sort,dist.between=dist.between,
                                  dist.part=dist.part,stable=stable,line=FALSE,arrow = FALSE,scale = FALSE,
                                  debug = FALSE,y=y.msa,arc.lty=conflict.lty,add = TRUE,col = conflict.col,x = x,lwd = conflict.lwd)
        }
      }else{
        warning("No conflicting interactions")
      }
    }

    if(nrow(helix.notconflict) != 0){
      helix.notconflict <- copy(helix.notconflict[,!"conflict"])
      helix.notconflict <- as.data.table(helix.notconflict)
      helix.notconflict <- expandHelixMltp(helix.notconflict)
      attr(helix.notconflict,"length") <- attr(helix,"length")
      plotHelixMltpDoubleLine(helix.notconflict,top.name = top.name,sort =sort,dist.between=dist.between,
                              dist.part=dist.part,stable=stable,line=FALSE,arrow = FALSE,scale = FALSE,
                              debug = FALSE,y=y.msa,arc.lty=1,add = TRUE,x = x)
    }else{
      warning("\nNo not conflict interactions.
              Use conflict.col and conflict.lty to display conflict interactions")
    }
  }


  #plot msa part
  for(i in 1:length(msa.cols.storage)){
    if(i ==1 ){
      if(is.vector((msa.seq.storage[[i]]))){
        if (grid) {
          plotCovarianceGrid(uniq[i,start] + 1, uniq[i,y] - height.msa * j, msa.cols.storage[[i]], height = height.msa,
                             border = grid.col, lwd = grid.lwd)
          if (text) {
            plotCovarianceText(uniq[i,start] + 1, uniq[i,y] - height.msa * j, msa.seq.storage[[i]], height = height.msa,
                               cex = text.cex, font = text.font, family = text.family,
                               col = text.col)
          }
          if (species > 0) {
            plotCovarianceSpecies(species - uniq[i,start] + 1, uniq[i,y] - height.msa * j,
                                  msa.name.storage[[i]], height = height.msa, cex = species.cex,
                                  font = species.font, family = species.family,
                                  col = species.col)
          }
        } else {
          plotCovarianceLine(uniq[i,start] + 1, uniq[i,y] - j, msa.cols.storage[[i]])
        }

      }else{
        for(j in 1:nrow(msa.seq.storage[[i]])) {
          if (grid) {
            plotCovarianceGrid(uniq[i,start] + 1, uniq[i,y] - height.msa * j, msa.cols.storage[[i]][j,], height = height.msa,
                               border = grid.col, lwd = grid.lwd)
            if (text) {
              plotCovarianceText(uniq[i,start] + 1, uniq[i,y] - height.msa * j, msa.seq.storage[[i]][j,], height = height.msa,
                                 cex = text.cex, font = text.font, family = text.family,
                                 col = text.col)
            }
            if (species > 0) {
              plotCovarianceSpecies(species - uniq[i,start] + 1, uniq[i,y] - height.msa * j,
                                    msa.name.storage[[i]][j], height = height.msa, cex = species.cex,
                                    font = species.font, family = species.family,
                                    col = species.col)
            }
          } else {
            plotCovarianceLine(uniq[i,start] + 1, uniq[i,y] - j, msa.cols.storage[[i]][j,])
          }
        }
      }

    }else{
      if(is.vector((msa.seq.storage[[i]]))){
        if (grid) {
          plotCovarianceGrid(uniq[i,start] + 1, uniq[i,y] + height.msa * j-1, msa.cols.storage[[i]], height = height.msa,
                             border = grid.col, lwd = grid.lwd)
          if (text) {
            plotCovarianceText(uniq[i,start] + 1, uniq[i,y] + height.msa * j-1, msa.seq.storage[[i]], height = height.msa,
                               cex = text.cex, font = text.font, family = text.family,
                               col = text.col)
          }
          if (species > 0) {
            plotCovarianceSpecies(species - uniq[i,start] + 1, uniq[i,y] + height.msa * j-1,
                                  msa.name.storage[[i]], height = height.msa, cex = species.cex,
                                  font = species.font, family = species.family,
                                  col = species.col)
          }
        } else {
          plotCovarianceLine(uniq[i,start] + 1, uniq[i,y] + j-1, msa.cols.storage[[i]])
        }

      }else{
        for(j in 1:nrow(msa.seq.storage[[i]])) {
          if (grid) {
            plotCovarianceGrid(uniq[i,start] + 1, uniq[i,y] + height.msa * j-1, msa.cols.storage[[i]][j,], height = height.msa,
                               border = grid.col, lwd = grid.lwd)
            if (text) {
              plotCovarianceText(uniq[i,start] + 1, uniq[i,y] + height.msa * j-1, msa.seq.storage[[i]][j,], height = height.msa,
                                 cex = text.cex, font = text.font, family = text.family,
                                 col = text.col)
            }
            if (species > 0) {
              plotCovarianceSpecies(species - uniq[i,start] + 1, uniq[i,y] + height.msa * j-1,
                                    msa.name.storage[[i]][j], height = height.msa, cex = species.cex,
                                    font = species.font, family = species.family,
                                    col = species.col)
            }
          } else {
            plotCovarianceLine(uniq[i,start] + 1, uniq[i,y] + j-1, msa.cols.storage[[i]][j,])
          }
        }

      }
    }
  }




  #plot legend
  if (legend) {
    if (base.colour) {
      text.legend <- c("A", "U", "G","C", "-", "?","T")
    } else {
      text.legend <- c("Conservation", "Covariation", "One-sided", "Invalid",
                       "Unpaired", "Gap", "Ambiguous")
    }
    legend("bottom", text.legend[used], border = NA, fill = palette[used],
           horiz = TRUE, bty = "n", xpd = NA)
  }
  #restore warnings
  options(warn = 1)

}


#' @rdname plotCovarianceMltpSingleLine
#' @export
#'
##########################################################
#plotComparisonHelixMltpDoubleLine
#########################################################
plotComparisonHelixMltpDoubleLine <- function(helix1,helix2,scale = TRUE, top.name = "FALSE", sort = TRUE,
                                              dist.between = 0.1, dist.part = 0.2,dist.par.plot=0.2, stable = FALSE, line = TRUE,
                                              arrow = TRUE, col.line = "black", col.arrow = "black",
                                              col = "black", shape = "circle",debug = FALSE, append = TRUE,
                                              add = FALSE, cols = NA, breaks = NA,log = FALSE,
                                              include.lowest = TRUE,...){


  options(warn = -1)
  if (!is.helix.mltp(helix1)) {
    stop("Input not a valid helix data frame, aborting")
  }
  if (!is.helix.mltp(helix2)) {
    stop("Input not a valid helix data frame, aborting")
  }
  if (unique(helix1$length == 1) != 1) {
    helix1 <- expandHelixMltp(helix1)
  }
  if (unique(helix2$length == 1) != 1) {
    helix2 <- expandHelixMltp(helix2)
  }


  dfTP <- merge(x = helix1, y = helix2, by = c("i", "j", "length", "ent.i", "ent.j"))

  if ("value.x" %in% names(dfTP)) {
    if ("value.y" %in% names(dfTP)) {
      dfTP <- dfTP[, 1:6]
      colnames(dfTP) <- c("i", "j", "length", "ent.i", "ent.j", "value")
    }
    else {
      colnames(dfTP) <- c("i", "j", "length","ent.i", "ent.j", "value")
    }
  }

  dfTP <- dfTP[, c("i", "j", "length", "value", "ent.i", "ent.j")]
  dfFN <- helix1[!dfTP, on = c("i", "j", "length","ent.i", "ent.j")]
  dfFP <- helix2[!dfTP, on = c("i", "j", "length", "ent.i", "ent.j")]

  if (nrow(dfTP) == 0) {
    stop("No True Positives")
  }
  if (nrow(dfFN) == 0) {
    stop("No False Negatives")
  }
  if (nrow(dfFP) == 0) {
    stop("No False Positives")
  }
  attr(dfTP, "length") <- attr(helix1, "length")
  attr(dfFN, "length") <- attr(helix1, "length")
  attr(dfFN, "length") <- attr(helix1, "length")

  if(unique(dfTP$value) ==1 ){
    warning("True Positives have same value\n Colour to Blue")
    dfTP[,col := "#4169E1"]
  }else{
    if(is.na(breaks)){
      if(is.na(cols)){
        dfTP <- colourByValueMltp(dfTP,get = TRUE, log = log, include.lowest = include.lowest)
      }else{
        dfTP <- colourByValueMltp(dfTP,get = TRUE, cols = cols, log = log, include.lowest = include.lowest)
      }
    }else{
      if(is.na(cols)){
        dfTP <- colourByValueMltp(dfTP,get = TRUE,breaks = breaks, log = log, include.lowest = include.lowest)
      }else{
        dfTP <- colourByValueMltp(dfTP,get = TRUE,breaks = breaks, cols = cols, log = log, include.lowest = include.lowest)
      }
    }

  }
  if(unique(dfFP$value) ==1 ){
    warning("False Positives have same value\n Colour to Blue")
    dfFP[,col := "#4169E1"]
  }else{
    if(is.na(breaks)){
      if(is.na(cols)){
        dfFP <- colourByValueMltp(dfFP,get = TRUE, log = log, include.lowest = include.lowest)
      }else{
        dfFP <- colourByValueMltp(dfFP,get = TRUE, cols = cols, log = log, include.lowest = include.lowest)
      }
    }else{
      if(is.na(cols)){
        dfFP <- colourByValueMltp(dfFP,get = TRUE,breaks = breaks, log = log, include.lowest = include.lowest)
      }else{
        dfFP <- colourByValueMltp(dfFP,get = TRUE,breaks = breaks, cols = cols, log = log, include.lowest = include.lowest)
      }
    }
  }


  #left side dfTP(coloured) and dfFN(full black)
  dfFN$col <- "#000000"
  df.left <- rbind(dfTP,dfFN[order(value)])
  attr(df.left, "length") <- attr(helix1, "length")


  #right side dfFP
  df.right <- dfFP
  attr(df.right, "length") <- attr(helix1, "length")


  height1 <- max.height.double(helix = df.left, comp.seq = comp.seq,top.name = top.name,
                               sort = sort, stable = stable,dist.between = dist.between,
                               dist.part = dist.part)

  height2 <- max.height.double(helix = df.right, comp.seq = comp.seq,top.name = top.name,
                               sort = sort, stable = stable,dist.between = dist.between,
                               dist.part = dist.part)
  width1 <- max.width.double(helix = df.left,comp.seq = comp.seq,top.name = top.name,sort = sort,
                             stable = stable,dist.between = dist.between, dist.part = dist.part)

  width2 <- max.width.double(helix = df.right,comp.seq = comp.seq,top.name = top.name,sort = sort,
                             stable = stable,dist.between = dist.between, dist.part = dist.part)


  blankPlotMod(width = width1+width2+(max(width1,width2)*dist.par.plot),
               top = max(height1[1],height2[1]),bottom = min(height1[2],height2[2]),
               scale = scale,debug = debug,...)


  plotHelixMltpDoubleLine(df.left, top.name = top.name, sort = sort,
                          dist.between = dist.between, dist.part = dist.part, stable = stable,
                          line = line, arrow = arrow, col.line = col.line, col.arrow = col.arrow,
                          col = col, shape = shape, scale = scale, debug = debug,
                          add = TRUE, append = FALSE,...)
  dist.part = dist.part
  dist.between = dist.between
  plotHelixMltpDoubleLine(df.right, top.name = top.name, sort = sort,
                          dist.between = dist.between, dist.part = dist.part, stable = stable,
                          line = line, arrow = arrow, col.line = col.line, col.arrow = col.arrow,
                          col = col, shape = shape, scale = scale, debug = debug,
                          add = TRUE, append = FALSE,x = width1 + (max(width1,width2)*dist.par.plot))


  options(warn = 1)

}



#' @rdname plotCovarianceMltpSingleLine
#' @export
#'
#########################################################
#plotCovarianceComparisonMltpDoubleLine
#########################################################
plotCovarianceComparisonMltpDoubleLine <- function(helix1,helix2,msa,
                                                   comp.seq = FALSE, y = 0, msa.all.seq = FALSE,
                                                   arcs = TRUE, top.name = "FALSE", sort = FALSE, conflict.cutoff = 0.01,
                                                   stable = FALSE, grid = TRUE, dist.between = 0.1, dist.part = 0.2,dist.par.plot=0.2,
                                                   add = FALSE, append = FALSE, pad = c(0, 0, 0, 0), shape = "circle",scale = FALSE,
                                                   debug = FALSE, col = "black",
                                                   base.colour = FALSE, palette = NA, text = TRUE, grid.col = "white",
                                                   grid.lwd = 0, text.cex = 0.2, text.col = "gray", text.font = 2,
                                                   text.family = "sans", species.cex = 0.2, species.col = "black",
                                                   species.font = 2, species.family = "mono", species = 10,
                                                   legend = TRUE, col.line = "black", col.arrow = "black",
                                                   cols = NA, breaks= NA,log = FALSE, include.lowest = TRUE,...){

  options(warn = -1)
  if (!is.helix.mltp(helix1)) {
    stop("Input not a valid helix data frame, aborting")
  }
  if (!is.helix.mltp(helix2)) {
    stop("Input not a valid helix data frame, aborting")
  }
  if (unique(helix1$length == 1) != 1) {
    helix1 <- expandHelixMltp(helix1)
  }
  if (unique(helix2$length == 1) != 1) {
    helix2 <- expandHelixMltp(helix2)
  }


  dfTP <- merge(x = helix1, y = helix2, by = c("i", "j", "length", "ent.i", "ent.j"))

  if ("value.x" %in% names(dfTP)) {
    if ("value.y" %in% names(dfTP)) {
      dfTP <- dfTP[, 1:6]
      colnames(dfTP) <- c("i", "j", "length", "ent.i", "ent.j", "value")
    }
    else {
      colnames(dfTP) <- c("i", "j", "length","ent.i", "ent.j", "value")
    }
  }

  dfTP <- dfTP[, c("i", "j", "length", "value", "ent.i", "ent.j")]
  dfFN <- helix1[!dfTP, on = c("i", "j", "length","ent.i", "ent.j")]
  dfFP <- helix2[!dfTP, on = c("i", "j", "length", "ent.i", "ent.j")]

  if (nrow(dfTP) == 0) {
    stop("No True Positives")
  }
  if (nrow(dfFN) == 0) {
    stop("No False Negatives")
  }
  if (nrow(dfFP) == 0) {
    stop("No False Positives")
  }
  attr(dfTP, "length") <- attr(helix1, "length")
  attr(dfFN, "length") <- attr(helix1, "length")
  attr(dfFN, "length") <- attr(helix1, "length")

  if(unique(dfTP$value) ==1 ){
    warning("True Positives have same value\n Colour to Blue")
    dfTP[,col := "#4169E1"]
  }else{
    if(is.na(breaks)){
      if(is.na(cols)){
        dfTP <- colourByValueMltp(dfTP,get = TRUE, log = log, include.lowest = include.lowest)
      }else{
        dfTP <- colourByValueMltp(dfTP,get = TRUE, cols = cols, log = log, include.lowest = include.lowest)
      }
    }else{
      if(is.na(cols)){
        dfTP <- colourByValueMltp(dfTP,get = TRUE,breaks = breaks, log = log, include.lowest = include.lowest)
      }else{
        dfTP <- colourByValueMltp(dfTP,get = TRUE,breaks = breaks, cols = cols, log = log, include.lowest = include.lowest)
      }
    }

  }
  if(unique(dfFP$value) ==1 ){
    warning("False Positives have same value\n Colour to Blue")
    dfFP[,col := "#4169E1"]
  }else{
    if(is.na(breaks)){
      if(is.na(cols)){
        dfFP <- colourByValueMltp(dfFP,get = TRUE, log = log, include.lowest = include.lowest)
      }else{
        dfFP <- colourByValueMltp(dfFP,get = TRUE, cols = cols, log = log, include.lowest = include.lowest)
      }
    }else{
      if(is.na(cols)){
        dfFP <- colourByValueMltp(dfFP,get = TRUE,breaks = breaks, log = log, include.lowest = include.lowest)
      }else{
        dfFP <- colourByValueMltp(dfFP,get = TRUE,breaks = breaks, cols = cols, log = log, include.lowest = include.lowest)
      }
    }
  }


  #left side dfTP(coloured) and dfFN(full black)
  dfFN$col <- "#000000"
  df.left <- rbind(dfTP,dfFN)
  attr(df.left, "length") <- attr(helix1, "length")


  #right side dfFP
  df.right <- dfFP
  attr(df.right, "length") <- attr(helix1, "length")


  height1 <- max.height.double(helix = df.left, msa = msa,comp.seq = comp.seq,
                               msa.all.seq = msa.all.seq,top.name = top.name,sort = sort,
                               conflict.cutoff = conflict.cutoff, stable = stable,
                               dist.between = dist.between, dist.part = dist.part,pad = pad,
                               species = species,arcs = arcs, palette = NA,base.colour = base.colour,
                               grid = grid)

  height2 <- max.height.double(helix = df.right, msa = msa,comp.seq = comp.seq,
                               msa.all.seq = msa.all.seq,top.name = top.name,sort = sort,
                               conflict.cutoff = conflict.cutoff, stable = stable,
                               dist.between = dist.between, dist.part = dist.part,pad = pad,
                               species = species,arcs = arcs, palette = NA,base.colour = base.colour,
                               grid = grid)

  width1 <- max.width.double(helix = df.left, msa = msa,comp.seq = comp.seq,
                             msa.all.seq = msa.all.seq,top.name = top.name,
                             sort = sort, conflict.cutoff = conflict.cutoff,
                             stable = stable,dist.between = dist.between,
                             dist.part = dist.part,pad = pad, species = species,
                             arcs = arcs, palette = NA)

  width2 <- max.width.double(helix = df.right, msa = msa,comp.seq = comp.seq,
                             msa.all.seq = msa.all.seq,top.name = top.name,
                             sort = sort, conflict.cutoff = conflict.cutoff,
                             stable = stable,dist.between = dist.between,
                             dist.part = dist.part,pad = pad, species = species,
                             arcs = arcs, palette = NA)


  blankPlotMod(width = width1+width2+(max(width1,width2)*dist.par.plot)+species,
               top = max(height1[1],height2[1]),bottom = min(height1[2],height2[2]),
               debug = debug,scale = scale,...)


  ################################################################

  #####Extract sequences and make msa as 1 BStringSet
  #####For next using on default functions
  #transfor fasta to data.table format
  msa.data.table <- bstring.to.data.table(msa,comp.seq = comp.seq)
  #call uniq from helix file
  uniq.fake <- copy(attr(helix1,"length"))

  #check if N of unique entities are same
  if(length(unique(uniq.fake[,group])) != length(unique(msa.data.table[,group]))){
    stop("Not equal amount of entities between fasta and msa")
  }
  if(!unique(sort(unique(msa.data.table[,group])) == sort(unique(msa.data.table[,group])))){
    stop("Entities names are not same between fasta and msa")
  }

  uniq.fake <- uniq.fake[,("width"):= lapply(.SD, as.numeric), .SDcols = "width"]
  index <- 1:uniq.fake[,.N]
  uniq.fake[,N:=index]
  #update uniq if it needed by top.name
  uniq.fake <- update.uniq.by.group.single(uniq.fake,top.name,sort)

  #contain full msa sequences
  msa.group.seq.full <- list()
  #contain names of sequences
  seq.names.msa <- list()
  for(i in 1:nrow(uniq.fake)){
    msa.group.seq.full[[uniq.fake[i,group]]] <- msa.data.table[group %like% uniq.fake[i,group]][,!"group"][,!"rest"]
    if(comp.seq){
      seq.names.msa[[uniq.fake[i,group]]] <- msa.group.seq.full[[uniq.fake[i,group]]][,seq.name]
    }
    if(!comp.seq){
      seq.names.msa[[uniq.fake[i,group]]] <- msa.group.seq.full[[uniq.fake[i,group]]][,name]
    }

  }
  #define duplicates and take massage to user, after delete duplicate names
  for(i in 1:nrow(uniq.fake)){
    if(length(unique(duplicated(seq.names.msa[[uniq.fake[i,group]]])))>1){
      message("\nDuplicate sequence name present.\nJust first sequence will be used")
    }
    seq.names.msa[[uniq.fake[i,group]]] <- unique(seq.names.msa[[uniq.fake[i,group]]])
  }

  if(comp.seq){
    #find common names between all entities
    common.seq.names <- Reduce(intersect, seq.names.msa)
    if(length(common.seq.names)==0){
      warning("No similar seq names present in entities\nFirst will be choosen")
    }
  }
  #define minimum of msa seq
  msa.rows.height <- c()
  for(i in 1:nrow(uniq.fake)){
    msa.rows.height <-c(msa.rows.height,nrow(msa.group.seq.full[[uniq.fake[i,group]]]))
  }
  height.msa.min <- min(msa.rows.height)
  height.msa.max <- max(msa.rows.height)

  #get height of msa
  height.msa <- ifelse(grid,1.25,1)
  #chose height by amount of plotted sequences
  if(msa.all.seq){
    msa.size <- height.msa*height.msa.max
  }else{
    msa.size <- height.msa*height.msa.min
  }


  plotHelixMltpDoubleLine(df.left, top.name = top.name, sort = sort, dist.between = dist.between,
                          dist.part = dist.part, stable = stable, line = FALSE, arrow = FALSE,
                          col.line = col.line, col.arrow = col.arrow, col = col,
                          shape = shape, scale = scale, debug = debug, add = TRUE,
                          y = y+msa.size+1, append = FALSE)

  plotHelixMltpDoubleLine(df.right, top.name = top.name, sort = sort, dist.between = dist.between,
                          dist.part = dist.part, stable = stable, line = FALSE, arrow = FALSE,
                          col.line = col.line, col.arrow = col.arrow, col = col,
                          shape = shape, scale = scale, debug = debug, add = TRUE,
                          y = y+msa.size+1, append = FALSE, x = width1 + (max(width1,width2)*dist.par.plot)+ species)

  ################################################################

  plotCovarianceMltpDoubleLine(helix = df.left, msa = msa,
                               comp.seq = comp.seq, msa.all.seq = msa.all.seq, arcs = FALSE,
                               top.name = top.name, sort = sort, conflict.cutoff = conflict.cutoff,
                               stable = stable, grid = grid, dist.between = dist.between,
                               dist.part = dist.part, add = TRUE, append = FALSE, pad = pad,
                               shape = shape, debug = debug, col = col, conflict.col = conflict.col,
                               conflict.lty = conflict.lty,conflict.lwd = conflict.lwd, base.colour = base.colour,
                               palette = palette, text = text, grid.col = grid.col,
                               grid.lwd = grid.lwd, text.cex = text.cex, text.col = text.col,
                               text.font = text.font, text.family = text.family, species.cex = species.cex,
                               species.col = species.col, species.font = species.font,
                               species.family = species.family, species = species, legend = FALSE,
                               col.line = col.line, col.arrow = col.arrow)

  plotCovarianceMltpDoubleLine(helix = df.right, msa = msa,
                               comp.seq = comp.seq, msa.all.seq = msa.all.seq, arcs = FALSE,
                               top.name = "FALSE", sort = sort, conflict.cutoff = conflict.cutoff,
                               stable = stable, grid = grid, dist.between = dist.between,
                               dist.part = dist.part, add = TRUE, append = FALSE, pad = pad,
                               shape = shape, debug = debug, col = col, conflict.col = conflict.col,
                               conflict.lty = conflict.lty, conflict.lwd = conflict.lwd, base.colour = base.colour,
                               palette = palette, text = text, grid.col = grid.col,
                               grid.lwd = grid.lwd, text.cex = text.cex, text.col = text.col,
                               text.font = text.font, text.family = text.family, species.cex = species.cex,
                               species.col = species.col, species.font = species.font,
                               species.family = species.family, species = species, legend = FALSE,
                               col.line = col.line, col.arrow = col.arrow,x = width1 + (max(width1,width2)*dist.par.plot)+ species)


  if (all(is.na(palette))) {
    if (base.colour) {
		  palette <- c("#E41A1C", "#4DAF4A", "#FF7F00", "#377EB8", "#BDBDBD","#984EA3","#00008b")
    }else {
      palette <- c("#00A651", "#0072BC", "#00B9F2",
                   "#F15A22", "#231F20", "#AAAAAA",
                   "#DA6FAB")
    }
  }else {
    if (length(palette) < ifelse(base.colour, 5, 7)) {
      palette <- rep(palette, length.out = ifelse(base.colour,5, 7))
    }
  }
  if (legend) {
    if (base.colour) {
      text.legend <- c("A", "U", "G","C", "-", "?","T")
    }else {
      text.legend <- c("Conservation", "Covariation",
                       "One-sided", "Invalid", "Unpaired",
                       "Gap", "Ambiguous")
    }
    legend("bottomleft", text.legend, border = NA,
           fill = palette, horiz = TRUE, bty = "n", xpd = NA)
  }


}


#' @rdname plotCovarianceMltpSingleLine
#' @export
#'
###############################################################################
#plotCovarianceDoubleMltpDoubleLine
###############################################################################
plotCovarianceDoubleMltpDoubleLine <- function(helix.left,helix.right,msa.left,msa.right,comp.seq = FALSE,
                                               msa.all.seq = TRUE,arcs = TRUE,
                                               top.name = "FALSE",sort = FALSE,scale = FALSE,
                                               conflict.cutoff = 0.01,stable = FALSE,grid = TRUE,
                                               dist.between = .1,dist.part = .2,dist.par.plot = .2,
                                               add = FALSE,append=FALSE,pad = c(0, 0, 0, 0),
                                               shape = "circle",debug = FALSE,
                                               col = "black",conflict.col = "blue",
                                               conflict.lty = 2,conflict.lwd = 1,base.colour = FALSE,
                                               palette =NA,text = TRUE,grid.col = "white",
                                               grid.lwd = 0,text.cex = 0.2,text.col = "gray",
                                               text.font = 2,text.family = "sans",
                                               species.cex = 0.2,species.col = "black",
                                               species.font = 2,species.family = "mono",
                                               species = 10,legend = TRUE,col.line = "black",
                                               col.arrow = "black",...){

  height1 <- max.height.double(helix = helix.left, msa = msa.left,comp.seq = comp.seq,
                               msa.all.seq = msa.all.seq,top.name = top.name,sort = sort,
                               conflict.cutoff = conflict.cutoff, stable = stable,
                               dist.between = dist.between, dist.part = dist.part,pad = pad,
                               species = species,arcs = arcs, palette = NA,base.colour = base.colour,
                               grid = grid)

  height2 <- max.height.double(helix = helix.right, msa = msa.right,comp.seq = comp.seq,
                               msa.all.seq = msa.all.seq,top.name = top.name,sort = sort,
                               conflict.cutoff = conflict.cutoff, stable = stable,
                               dist.between = dist.between, dist.part = dist.part,pad = pad,
                               species = species,arcs = arcs, palette = NA,base.colour = base.colour,
                               grid = grid)

  width1 <- max.width.double(helix = helix.left, msa = msa.left,comp.seq = comp.seq,
                             msa.all.seq = msa.all.seq,top.name = top.name,
                             sort = sort, conflict.cutoff = conflict.cutoff,
                             stable = stable,dist.between = dist.between,
                             dist.part = dist.part,pad = pad, species = species,
                             arcs = arcs, palette = NA)

  width2 <- max.width.double(helix = helix.right, msa = msa.right,comp.seq = comp.seq,
                             msa.all.seq = msa.all.seq,top.name = top.name,
                             sort = sort, conflict.cutoff = conflict.cutoff,
                             stable = stable,dist.between = dist.between,
                             dist.part = dist.part,pad = pad, species = species,
                             arcs = arcs, palette = NA)


  blankPlotMod(width = width1+width2+(max(width1,width2)*dist.par.plot)+species,
               top = max(height1[1],height2[1]),bottom = min(height1[2],height2[2]),
               scale = scale,debug = debug,...)



  plotCovarianceMltpDoubleLine(helix = helix.left,msa = msa.left, comp.seq = comp.seq,
                               msa.all.seq = msa.all.seq,arcs = arcs,
                               top.name = top.name,sort = sort,
                               conflict.cutoff = conflict.cutoff,stable = stable,grid = grid,
                               dist.between = dist.between,dist.part = dist.part,
                               add = TRUE,append = FALSE,pad = pad,
                               shape = shape,debug = debug,
                               col = col,conflict.col = conflict.col,
                               conflict.lty = conflict.lty,conflict.lwd = conflict.lwd,base.colour = base.colour,
                               palette =palette,text = text,grid.col = grid.col,
                               grid.lwd = grid.lwd,text.cex = text.cex,text.col = text.col,
                               text.font = text.font,text.family = text.family,
                               species.cex = species.cex,species.col = species.col,
                               species.font = species.font,species.family = species.family,
                               species = species,legend = FALSE,col.line = col.line,
                               col.arrow = col.arrow,...)


  plotCovarianceMltpDoubleLine(helix = helix.right,msa = msa.right, comp.seq = comp.seq,
                               msa.all.seq = msa.all.seq,arcs = arcs,
                               top.name = "FALSE",sort = sort,
                               conflict.cutoff = conflict.cutoff,stable = stable,grid = grid,
                               dist.between = dist.between,dist.part = dist.part,
                               add = TRUE,append = FALSE,pad = pad,
                               shape = shape,debug = debug,
                               col = col,conflict.col = conflict.col,
                               conflict.lty = conflict.lty,conflict.lwd = conflict.lwd,base.colour = base.colour,
                               palette =palette,text = text,grid.col = grid.col,
                               grid.lwd = grid.lwd,text.cex = text.cex,text.col = text.col,
                               text.font = text.font,text.family = text.family,
                               species.cex = species.cex,species.col = species.col,
                               species.font = species.font,species.family = species.family,
                               species = species,legend = FALSE,col.line = col.line,
                               col.arrow = col.arrow,x = width1 + (max(width1,width2)*dist.par.plot)+ species)


  #pallete colours
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

  #plot legend
  if (legend) {
    if (base.colour) {
      text.legend <- c("A", "U", "G", "C", "-", "?")
    } else {
      text.legend <- c("Conservation", "Covariation", "One-sided", "Invalid",
                       "Unpaired", "Gap", "Ambiguous")
    }
    legend("bottomleft", text.legend, border = NA, fill = palette,
           horiz = TRUE, bty = "n", xpd = NA)
  }


}

