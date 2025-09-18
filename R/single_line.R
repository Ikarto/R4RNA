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


#' @export
#'
###############################################################################
#plotComparisonHelix
###############################################################################
plotComparisonHelix <- function(predict, known, miss = "black", line = TRUE,
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
    blankPlotMod(width, top, -bot, ...)
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






#' Plot Helix in arc diagram in single or double line mode
#'
#' @description
#'
#' Plots a helix data table as an arc diagram in single line or double line mode,
#' with styling possible with properly named additional columns on the data table.
#'
#' @param helix Helix data.tables, with the 6 mandatory columns. "col" refer to a styling column, and will be used for styling the helix. See example for styling usage.
#' @param top.name Name of entity, which will be set as top name and will refer to first entity in plot
#' @param sort sort entities by alphabetic order
#' @param dist.part distance between entities, dist.part*average length of entities
#' @param dist.between distance between entities in y axes, dist.between*average length of entities
#' @param stable used to have same height of resulted plot, based on length of entities
#' @param line If TRUE, a horizontal line representing the sequence is plotted.
#' @param arrow If TRUE, an arrow is played on the right end of the line.
#' @param col.line Colour of line, default is "black"
#' @param col.arrow Colour of arrow, default is "black"
#' @param col colour of arcs, if there is no style column in helix data.table sstructure
#' @param shape One of "circle", "triangle", or "square", specifying the shape of the arcs.
#' @param scale If TRUE, inserts a scale on the plot.
#' @param debug If TRUE, frames the boundaries of the intended plotting space in red, used to determine if inputs produce expected output area. Also outputs to STDIN dimensions of the plot.
#' @param add,append If TRUE, graphical elements are added to the active plot device, else a new plot device is created for the plot.
#' @param flip If TRUE, flips the arcs upside down about the y-axis.
#' @param y The vertical offset of the arc base relative to 0 along the y-axis.
#' @param arc.lty define lty of arcs
#'
#' @author Volodymyr Tsybulskyi
#' @seealso colourByCountMltp, colourByValueMltp
#' @details
#' plotHelixMltpSingleLine creates an arc diagram with all arcs on top and in single line mode,
#' entites goes from left to right.
#' plotHelixMltpDoubleLine creates an arc diagram with first entity below and other in top,
#' if you choose sort option then they will be sorted alphabetically. With top name you can define
#' which entity will in below.
#'
#' @export
#'
###############################################################################
#plotHelixMltpSingle
###############################################################################
plotHelixMltpSingleLine <- function(helix,top.name = "FALSE",sort = TRUE,
                                    dist.part = .2,stable = FALSE,line = TRUE,
                                    arrow = TRUE,col.line = "black",
                                    col.arrow = "black",col = "black",
                                    shape = "circle",scale = FALSE,
                                    debug = FALSE,add = FALSE,flip = FALSE,y=0,arc.lty = 1,...){
  #reduce data.table warning
  options(warn = -1)

  if(helix[length>1,.N]!=0){
    stop("Helix is not expanded")
  }
  if(is.null(attr(helix,"length"))){
    stop("No attribute for helix")
  }
  #read len data
  uniq <- attr(helix,"length")
  #sort and define top name
  uniq <- update.uniq.by.group.single(uniq,top.name,sort)
  #calculate values for uniq
  uniq <- prep.uniq.single(uniq,dist.part)
  full.len <- attr(uniq,"full.length")
  dist.between <- attr(uniq,"dist.between")
  #update helix to new coordinates
  helix <- update.helix.coord.single(helix,uniq)
  #define max height for not stable variant(dynamic height)
  max.height <- max.height.single(helix,full.len,stable)
  #create space
  if(!add){
    if(flip){
      blankPlotMod(full.len+dist.between/8, 1+y,-max.height+y,scale = scale,debug = debug,...)
    }else{
      blankPlotMod(full.len+dist.between/8, max.height+y, -1+y,scale = scale,debug = debug,...)
    }

  }
  y.uniq <- y
  uniq <- uniq[,y:=y.uniq]
  #plot bases with arrows if specified
  if(line){
    plot.bases.single(uniq,col.line,y=y)
  }
  if(arrow){
    plot.arrows.single(uniq,dist.between,col.arrow,y=y)
  }
  #plot interactions arcs
  plot.arcs.single(helix,col = col,flip = flip,y=y,lty=arc.lty,shape=shape)
  #restore warnings
  options(warn = 1)
}


#' Plot Helix in arc diagram in single line mode with double helix
#'
#' @description
#'
#' Plots a helix data table as an arc diagram in single line, where 2 helix are given by user.
#' with styling possible with properly named additional columns on the data table "col",
#' can be used separately in each helix data.table.
#'
#' @param top,bottom helixs data.table, with the six mandatory columns. "col" column can be specified as styling column.
#' @param top.name Name of entity, which will be set as top name and will refer to first entity in plot
#' @param sort sort entities by alphabetic order
#' @param dist.part distance between entities, dist.part*average length of entities
#' @param stable used to have same height of resulted plot, based on length of entities
#' @param line If TRUE, a horizontal line representing the sequence is plotted.
#' @param arrow If TRUE, an arrow is played on the right end of the line.
#' @param col.line Colour of line, default is "black"
#' @param col.arrow Colour of arrow, default is "black"
#' @param col colour of arcs, if there is no style column in helix data.table sstructure
#' @param shape One of "circle", "triangle", or "square", specifying the shape of the arcs.
#' @param scale If TRUE, inserts a scale on the plot.
#' @param debug If TRUE, frames the boundaries of the intended plotting space in red, used to determine if inputs produce expected output area. Also outputs to STDIN dimensions of the plot.
#' @param add If TRUE, graphical elements are added to the active plot device, else a new plot device is created for the plot.
#' @param y The vertical offset of the arc base relative to 0 along the y-axis.
#'
#' @author Volodymyr Tsybulskyi
#'
#' @export
#'
###############################################################################
#plotDoubleHelixMltpSingleLine
###############################################################################
plotDoubleHelixMltpSingleLine <- function(top,bottom,top.name = "FALSE",sort = TRUE,
                                          dist.part = .2,stable = FALSE,line = TRUE,
                                          arrow = TRUE,col.line = "black",
                                          col.arrow = "black",col = "black",
                                          shape = "circle",scale = FALSE,
                                          debug = FALSE,add = FALSE,y=0,...){

  #reduce data.table warning
  options(warn = -1)

  if(top[length>1,.N]!=0 | bottom[length>1,.N]!=0){
    stop("Helix is not expanded")
  }
  if(is.null(attr(top,"length")) | is.null(attr(bottom,"length"))){
    stop("No attribute for helix")
  }
  #check if length is same for both helices
  if(sum(attr(top,"length")[,as.numeric(width)]) !=
     sum(attr(bottom,"length")[,as.numeric(width)])){
    stop("Helices has no same length")
  }
  #######################################################
  #save colour if it present                            #
  #######################################################
  if("col" %in% names(top)){
    col.top <- top$col
  }
  if("col" %in% names(bottom)){
    col.bottom <- bottom$col
  }
  #######################################################
  #define width for plot,assume that both has same width#
  #######################################################
  #read len data
  uniq_1 <- attr(top,"length")
  uniq_2 <- attr(bottom,"length")
  #sort and define top name
  uniq_1 <- update.uniq.by.group.single(uniq_1,top.name,sort)
  uniq_2 <- update.uniq.by.group.single(uniq_2,top.name,sort)
  #calculate values for uniq
  uniq_1 <- prep.uniq.single(uniq_1,dist.part)
  uniq_2 <- prep.uniq.single(uniq_2,dist.part)
  full.len <- attr(uniq_1,"full.length")
  dist.between <- attr(uniq_1,"dist.between")
  #######################################################
  #define height for plot pos and neg                   #
  #######################################################
  #update helix to new coordinates
  helix_1 <- update.helix.coord.single(top,uniq_1)
  helix_2 <- update.helix.coord.single(bottom,uniq_2)

  max.height.top <- max.height.single(helix_1,full.len,stable)
  max.height.bottom <- max.height.single(helix_2,full.len,stable)
  #######################################################
  #plot data                                            #
  #######################################################
  if(!add){
    blankPlotMod(width = full.len+dist.between/8, top = max.height.top+y, bottom = -max.height.bottom+y,
                 scale = scale,debug = debug,...)
  }
  #######################################################
  #return col                                           #
  #######################################################
  if("col" %in% names(top)){
    top <- top[,1:6]
    top$col <- col.top
  }else{
    top <- top[,1:6]
  }

  if("col" %in% names(bottom)){
    bottom <- bottom[,1:6]
    bottom$col <- col.bottom
  }else{
    bottom <- bottom[,1:6]
  }

  plotHelixMltpSingleLine(top,col = col,flip = FALSE,line = FALSE,arrow = FALSE,add = TRUE,
                          scale = scale,dist.part = dist.part,stable = stable,
                          shape = shape,debug = debug,top.name = top.name,
                          sort = sort,y = y)
  plotHelixMltpSingleLine(bottom,col = col,flip = TRUE,line = FALSE,arrow = FALSE,add = TRUE,
                          scale = scale,dist.part = dist.part,stable = stable,
                          shape = shape,debug = debug,top.name = top.name,
                          sort = sort,y = y)

  uniq_1$y <- y

  if(line){
    plot.bases.single(uniq_1,col.line)
  }
  if(arrow){
    plot.arrows.single(uniq_1,dist.between,col.arrow)
  }
  #restore warnings
  options(warn = 1)
}




#'
#'Plot nucleotide sequence coloured by covariance
#'
#'@description
#'Given a multiple sequence alignment and a corresponding secondary structure,
#'nucleotides in the sequence alignment will be coloured according to the
#'basepairing and conservation status, where green is the most commonly
#'observed valid basepair in the column, dark blue being valid covariation
#'(i.e. mutation into another valid basepair), cyan is one-sided mutation that
#'retains the basepair, and red is a mutation where the basepair has been lost.
#'
#' @name plotCovarianceMltpSingleLine
#'
#'
#' @param msa,msa1,msa2,msa.left,msa.right Multiple sequence alignment as an array of named characters,
#' all of equal length. Typically output of readBstringSet. \code{msa.left} and \code{msa.right}
#' are specific to \code{helix.left} and \code{helix.right}. \code{msa1} and \code{msa2}
#' are specific to \code{helix1} and \code{helix2}
#' @param helix,helix1,helix2,helix.left,helix.right A helix data.table with a structure corresponding to entities in msa
#' @param comp.seq Option to compare sequences from different entities by their "name" and provide
#' proper covariance colouring in future
#' @param top.name Name of entity, which will be set as top name and will refer to first entity in plot
#' @param sort sort entities by alphabetic order
#' @param msa.all.seq plot all sequences, if number of
#' @param arcs TRUE if the structure should be plotted as arcs. Arcs may be
#' styled with styling columns, see example and plotHelixMltp for details.
#' @param stable used to have same height of resulted plot, based on length of entities
#' @param grid TRUE if the multiple sequence alignment is to be drawn as a
#' grid of bases, else the multiple sequence alignment is drawn as equidistant horizontal lines.
#' @param dist.part distance between entities, dist.part*average length of entities
#' @param dist.between distance between entities on y axes, dist.between*average length of entities
#' @param append If TRUE, graphical elements are added to the active plot device,
#' else a new plot device is created for the plot.
#' @param append If TRUE, graphical elements are added to the active plot device,
#' else a new plot device is created for the plot.
#' @param pad A four integer array passed to blankPlot, specifies the
#' number of pixels to pad the bottom, left, top and right sides of the figure with, repsectively.
#' @param shape One of "circle", "triangle", or "square", specifying the shape of the arcs.
#' @param scale If TRUE, inserts a scale on the plot.
#' @param debug If TRUE, frames the boundaries of the intended plotting space
#' in red, used to determine if inputs produce expected output area.
#' Also outputs to STDIN dimensions of the plot.
#' @param col colour of arcs, if there is no style column in helix data.table sstructure
#' @param conflict.col,conflict.lty Determines the line type (style)
#' and colour to be used for conflicting basepairs. By default,
#' conflicting helices are drawn as dotted lines (lty = 2) and
#' whatever colour was originally assigned to it (col = NA).
#' Conflicting helices may be coloured by setting conflict.col to
#' some R-compatible colour name. If both arguments are set to NA,
#' then no attempt to exclude conflicting helices will be made when
#' colouring covariance plot columns, which in most cases will
#' render the plot nonsensical.
#' @param base.colour TRUE if bases are to be coloured by
#' nucleotide instead of basepair conservation.
#' @param palette A list of colour names to override the
#' default colour palette. When base.colour is TRUE, the first
#' 6 colours will be used for colouring bases A, U, G, C, - (gap),
#' and ? (everything else), respectively. When base.colour is FALSE,
#' the first 7 colours will be used for colouring conserved basepairs,
#' covarying basepairs, one-sided conserved basepairs, invalid basepairs,
#' unpaired bases, gaps, and bases/pairs with ambiguous bases, resepctively.
#' If the palette is shorter than the expected length, the palette will
#' simply cycle. "NA" is a valid colour, that will effectively plot nothing.
#' @param text Only applicable when grid is TRUE.
#' TRUE if the grid is to be filled with nucleotide character.
#' @param grid.col,grid.lwd The colour and line width of the borders
#' displayed when grid is TRUE.
#' @param text.cex,text.col,text.font,text.family cex, col,
#' family and font for the text displayed via the text option. Use help("par")
#' for more information the paramters.
#' @param species.cex,species.col,species.font,species.family cex, col,
#' family and font for the species text displayed via the species option.
#' Use help("par") for more information the paramters.
#' @param species If a number greater than 0 is given,
#' then species names for the multiple sequence alignment
#' will be printed along the left side. This name is
#' typically the entire header lines of FASTA entries
#' from readFasta, and can be manually manipulated using
#' the names function. The number specifies the start
#' position relative to the left edge of the multiple sequence alignment).
#' @param legend TRUE if legend are to be shown.
#' @param col.line Colour of line, default is "black"
#' @param col.arrow Colour of arrow, default is "black"
#' @param y The vertical offset of the entire figure relative to 0.
#' @param dist.plot distance in pixels between two single line coavriances
#'
#' @export
#'
##############################################################
#plotCovarianceMltpSingleLine
##############################################################
plotCovarianceMltpSingleLine <- function(msa,helix,comp.seq = FALSE,top.name = "FALSE",sort = FALSE,
                                         arcs = TRUE,msa.all.seq = TRUE,conflict.cutoff = 0.01,stable = FALSE,
                                         y=0,flip = FALSE,grid = TRUE,dist.between = .1,
                                         dist.part = .2,add = FALSE,pad = c(0, 0, 0, 0),
                                         shape = "circle",debug = FALSE,scale = FALSE,
                                         col = "black",conflict.col = NA,conflict.lty = 2, conflict.lwd = 1,
                                         base.colour = FALSE,palette =NA,text = TRUE,
                                         grid.col = "white",grid.lwd = 0,text.cex = 0.2,
                                         text.col = "gray",text.font = 2,text.family = "sans",
                                         species.cex = 0.2,species.col = "black",
                                         species.font = 2,species.family = "mono",
                                         species = 10,legend = FALSE,...){

  if(!is.helix.mltp(helix)){
    stop("invalid helix input")
  }
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }
  #restore warnings
  options(warn = -1)

  shape.conflict <- shape

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


  #create blank plot
  if(arcs){
    max.height <- def.max.height.helix(helix,top.name = top.name,
                                       sort = sort,dist.part = dist.part,
                                       stable = stable)
    full.length <- attr(max.height,"full.length")
  }else{
    max.height <- def.max.height.helix(helix,top.name = top.name,
                                       sort = sort,dist.part = dist.part,
                                       stable = stable)
    full.length <- attr(max.height,"full.length")
    max.height <- 1
  }

  #get height of msa
  height.msa <- ifelse(grid,1.25,1)
  #chose height by amount of plotted sequences
  if(msa.all.seq){
    msa.size <- height.msa*height.msa.max
  }else{
    msa.size <- height.msa*height.msa.min
  }
  #define max height and
  if(flip){
    max.height.pos <- 1+y
    max.height.neg <-  -max.height-msa.size+y
  }else{
    max.height.pos <- max.height+y
    max.height.neg <- -msa.size+y
  }

  #plot blank plot
  if(!add){
    if(species > 0){
      pad[2] <- pad[2] + species
    }
    blankPlotMod(width = full.length, top = max.height.pos,bottom = max.height.neg, pad,debug = debug,scale = scale,...)
  }

  #update y coordinate for msa if flip option set upped
  if(flip){
    y <- y - msa.size
  }

  #plot conflict
  if(arcs){
    if (!is.na(conflict.col) || !is.na(conflict.lty)) {
      if(nrow(helix.conflict) != 0){

        helix.conflict <- copy(helix.conflict[,!"conflict"])
        helix.conflict <- as.data.table(helix.conflict)
        helix.conflict <- expandHelixMltp(helix.conflict)
        attr(helix.conflict,"length") <- attr(helix,"length")
        plotHelixMltpSingleLine(helix.conflict[,1:6],col = conflict.col,flip = flip,line = FALSE,arrow = FALSE,add = TRUE,
                                scale = FALSE,dist.part = dist.part,stable = stable,
                                shape = shape.conflict,debug = debug,top.name = top.name,
                                sort = sort,y = y,arc.lty=conflict.lty,lwd = conflict.lwd)
      }else{
        warning("No conflicting interactions")
      }
    }

    if(nrow(helix.notconflict) != 0){
      helix.notconflict <- copy(helix.notconflict[,!"conflict"])
      helix.notconflict <- as.data.table(helix.notconflict)
      helix.notconflict <- expandHelixMltp(helix.notconflict)
      attr(helix.notconflict,"length") <- attr(helix,"length")
      plotHelixMltpSingleLine(helix.notconflict,col = col,flip = flip,line = FALSE,arrow = FALSE,add = TRUE,
                              scale = FALSE,dist.part = dist.part,stable = stable,
                              shape = shape,debug = debug,top.name = top.name,
                              sort = sort,y = y)
    }else{
      warning("\nNo not conflict interactions.
              Use conflict.col and conflict.lty to display conflict interactions")
    }
  }

  #get start and end positions
  options(warn = -1)
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

  #extract separate group colors segments and create separate
  uniq.fake[,start:= start+1]
  #list with msa's
  msa.cols.storage <- list()
  msa.name.storage <- list()
  msa.seq.storage <- list()

  #msa.selected.list msa.group.seq.full
  for (i in 1:uniq.fake[,.N]) {
    if (is.vector(cols[,uniq.fake[i,start]:uniq.fake[i,end]])) {
      msa.cols.storage[[i]] <- matrix(cols[,uniq.fake[i,start]:uniq.fake[i,end]],nrow = 1)
    }else{
      msa.cols.storage[[i]] <- cols[,uniq.fake[i,start]:uniq.fake[i,end]]
    }

    if(nrow(msa.group.seq.full[[uniq.fake[i,group]]])!=nrow(msa.cols.storage[[i]])){
      for (j in 1:(nrow(msa.group.seq.full[[uniq.fake[i,group]]])-nrow(msa.cols.storage[[i]]))) {
        msa.cols.storage[[i]] <- rbind(msa.cols.storage[[i]],rep("#231F20",ncol(msa.cols.storage[[i]])))
      }
    }
    if(comp.seq){
      tmp <- rbind(msa.selected.list[[uniq.fake[i,group]]],msa.group.seq.full[[uniq.fake[i,group]]][!msa.selected.list[[uniq.fake[i,group]]],on = c("seq.name","sequence")])
    }else{
      tmp <- rbind(msa.selected.list[[uniq.fake[i,group]]],msa.group.seq.full[[uniq.fake[i,group]]][!msa.selected.list[[uniq.fake[i,group]]],on = c("name","sequence")])
    }

    tmp1 <- strsplit(tmp[,sequence],"")
    matrix <- matrix(ncol = length(tmp1[[1]]),nrow = length(tmp1))
    for (j in 1:length(tmp1)) {
      matrix[j,] <- tmp1[[j]]
    }
    msa.seq.storage[[i]] <- matrix
    if(comp.seq){
      msa.name.storage[[i]] <- tmp[,seq.name]
    }else{
      msa.name.storage[[i]] <- tmp[,name]
    }

  }


  #cut msa if it cutted version
  if(!msa.all.seq){
    for (i in 1:length(msa.cols.storage)) {
      msa.seq.storage[[i]] <- msa.seq.storage[[i]][1:height.msa.min,]
      msa.cols.storage[[i]] <- msa.cols.storage[[i]][1:height.msa.min,]
      msa.name.storage[[i]] <- msa.name.storage[[i]][1:height.msa.min]
    }
  }


  #plot msa part
  for(i in 1:length(msa.cols.storage)){
    for(j in 1:nrow(msa.seq.storage[[i]])) {
      if(flip){
        if (grid) {
          plotCovarianceGrid(uniq[i,start] + 1, y + height.msa * j-1, msa.cols.storage[[i]][j,], height = height.msa,
                             border = grid.col, lwd = grid.lwd)
          if (text) {
            plotCovarianceText(uniq[i,start] + 1, y + height.msa * j-1, msa.seq.storage[[i]][j,], height = height.msa,
                               cex = text.cex, font = text.font, family = text.family,
                               col = text.col)
          }
          if (species > 0) {
            plotCovarianceSpecies(species - uniq[i,start] + 1, y + height.msa * j-1,
                                  msa.name.storage[[i]][j], height = height.msa, cex = species.cex,
                                  font = species.font, family = species.family,
                                  col = species.col)
          }
        } else {
          plotCovarianceLine(uniq[i,start] + 1, y + j, msa.cols.storage[[i]][j,])
        }
      }else{
        if (grid) {
          plotCovarianceGrid(uniq[i,start] + 1, y - height.msa * j, msa.cols.storage[[i]][j,], height = height.msa,
                             border = grid.col, lwd = grid.lwd)
          if (text) {
            plotCovarianceText(uniq[i,start] + 1, y - height.msa * j, msa.seq.storage[[i]][j,], height = height.msa,
                               cex = text.cex, font = text.font, family = text.family,
                               col = text.col)
          }
          if (species > 0) {
            plotCovarianceSpecies(species - uniq[i,start] + 1, y - height.msa * j,
                                  msa.name.storage[[i]][j], height = height.msa, cex = species.cex,
                                  font = species.font, family = species.family,
                                  col = species.col)
          }
        } else {

          plotCovarianceLine(uniq[i,start] + 1, y - j, msa.cols.storage[[i]][j,])

        }
      }

    }
  }

  #plot legend
  if (legend) {
    if (base.colour) {
      text.legend <- c("A", "U", "G", "C", "-", "?","T")
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


#' @rdname plotHelixMltpSingleLine
#' @name plotComparisonHelixMltpSingleLine
#' @export
#'
###############################################################################
#plotComparisonHelixMltp
###############################################################################
plotComparisonHelixMltpSingleLine <- function(helix1,helix2,top.name = "FALSE",sort = TRUE,
                                              dist.between = .1,dist.part = .2,
                                              stable = FALSE,line = TRUE,arrow = TRUE,
                                              col.line = "black",col.arrow = "black",
                                              col = "black",shape = "circle",scale = FALSE,
                                              debug = FALSE,append=TRUE,add = FALSE,cols = NA,
                                              breaks = NA, log = FALSE, include.lowest = TRUE,...){

  #reduce data.table warning
  options(warn = -1)

  if (!is.helix.mltp(helix1)) {
    stop("Input not a valid helix data frame, aborting")
  }
  if (!is.helix.mltp(helix2)) {
    stop("Input not a valid helix data frame, aborting")
  }
  if(unique(helix1$length == 1) != 1){
    helix1 <- expandHelixMltp(helix1)
  }
  if(unique(helix2$length == 1) != 1){
    helix2 <- expandHelixMltp(helix2)
  }

  #True Positive (in prediction 1 and 2)
  dfTP <- merge(x=helix1,y=helix2,by=c("i","j","length","ent.i","ent.j" ))
  #update df
  if("value.x" %in% names(dfTP)){
    if("value.y" %in% names(dfTP)){
      dfTP <- dfTP[,1:6]
      colnames(dfTP) <- c("i","j","length","ent.i","ent.j","value")
    }else{
      colnames(dfTP) <- c("i","j","length","ent.i","ent.j","value")
    }
  }

  dfTP <- dfTP[,c("i","j","length","value","ent.i","ent.j")]

  #False Negative (just in prediction 1)
  dfFN <- helix1[!dfTP, on= c("i","j","length","ent.i","ent.j" )]
  #False Positive (just in prediction 2)
  dfFP <- helix2[!dfTP, on= c("i","j","length","ent.i","ent.j" )]

  if(nrow(dfTP) == 0){
    warning("No True Positives")
  }else{
    attr(dfTP, "length") <- attr(helix1,"length")

    if(length(unique(dfTP$value)) == 1){
      warning("True Positives have same value\n Colour to Blue")
      dfTP[,col := "#4169E1"]
    }else{
      if(is.na(breaks)){
        # if(is.na(cols)){
        if (missing(cols) || (length(cols) == 1 && is.na(cols))) {
          dfTP <- colourByValueMltp(dfTP,get = TRUE, log = log, include.lowest = include.lowest)
        }else{
          dfTP <- colourByValueMltp(dfTP,get = TRUE, cols = cols, log = log, include.lowest = include.lowest)
        }
      }else{
        # if(is.na(cols)){
        if (missing(cols) || (length(cols) == 1 && is.na(cols))) {
          dfTP <- colourByValueMltp(dfTP,get = TRUE,breaks = breaks, log = log, include.lowest = include.lowest)
        }else{
          dfTP <- colourByValueMltp(dfTP,get = TRUE,breaks = breaks, cols = cols, log = log, include.lowest = include.lowest)
        }
      }

    }
  }
  if(nrow(dfFN) == 0){
    warning("No False Negatives")
  }else{
    attr(dfFN, "length") <- attr(helix1,"length")

    if(length(unique(dfFP$value)) ==1 ){
      warning("False Positives have same value\n Colour to Blue")
      dfFP[,col := "#4169E1"]
    }else{
      if(is.na(breaks)){
        # if(is.na(cols)){
        if (missing(cols) || (length(cols) == 1 && is.na(cols))) {
          dfFP <- colourByValueMltp(dfFP,get = TRUE, log = log, include.lowest = include.lowest)
        }else{
          dfFP <- colourByValueMltp(dfFP,get = TRUE, cols = cols, log = log, include.lowest = include.lowest)
        }
      }else{
        # if(is.na(cols)){
        if (missing(cols) || (length(cols) == 1 && is.na(cols))) {
          dfFP <- colourByValueMltp(dfFP,get = TRUE,breaks = breaks, log = log, include.lowest = include.lowest)
        }else{
          dfFP <- colourByValueMltp(dfFP,get = TRUE,breaks = breaks, cols = cols, log = log, include.lowest = include.lowest)
        }
      }
    }
  }

  if(nrow(dfFP) == 0){
    warning("No False Positives")
  }else{
    attr(dfFP, "length") <- attr(helix1,"length")
  }

  if(nrow(dfTP) != 0){
    dfBP <- rbind(dfTP[,1:6])
    attr(dfBP, "length") <- attr(helix1,"length")
  }

  if(nrow(dfFN) != 0){
    dfBP <- rbind(dfFN[,1:6])
    attr(dfBP, "length") <- attr(helix1,"length")
  }

  if(nrow(dfFN) != 0 & nrow(dfTP) != 0){
    dfBP <- rbind(dfFN[,1:6],dfTP[,1:6])
    attr(dfBP, "length") <- attr(helix1,"length")
  }

  max.height.pos <- def.max.height.helix(dfBP,top.name = top.name,
                                         sort = sort,dist.part = dist.part,
                                         stable = stable)
  if(nrow(dfFP)==0){
    max.height.neg <- -max.height.pos
  }else{
    max.height.neg <- -def.max.height.helix(dfFP,top.name = top.name,
                                            sort = sort,dist.part = dist.part,
                                            stable = stable)
  }

  dist.between <- attr(max.height.pos,"dist.between")
  uniq <- attr(max.height.pos,"uniq")
  full.len <- uniq[.N,end]

  if(!add){
    blankPlotMod(full.len+dist.between/8, max.height.pos, max.height.neg,
                 scale = scale,debug = debug,...)
  }

  if(nrow(dfFN)!=0){
    plotHelixMltpSingleLine(dfFN,col = col,flip = FALSE,line = FALSE,arrow = FALSE,add = TRUE,
                            scale = scale,dist.part = dist.part,stable = stable,
                            shape = shape,debug = debug,top.name = top.name,
                            sort = sort,...)
  }
  if(nrow(dfTP)!=0){
    plotHelixMltpSingleLine(dfTP[order(value)],col = col,flip = FALSE,line = FALSE,arrow = FALSE,add = TRUE,
                            scale = scale,dist.part = dist.part,stable = stable,
                            shape = shape,debug = debug,top.name = top.name,
                            sort = sort,...)
  }
  if(nrow(dfFP)!=0){
    plotHelixMltpSingleLine(dfFP[order(value)],col = col,flip = TRUE,line = FALSE,arrow = FALSE,add = TRUE,
                            scale = scale,dist.part = dist.part,stable = stable,
                            shape = shape,debug = debug,top.name = top.name,
                            sort = sort,...)
  }

  uniq$y <- 0

  if(line){
    plot.bases.single(uniq,col.line)
  }
  if(arrow){
    plot.arrows.single(uniq,dist.between,col.arrow)
  }

  #restore warnings
  options(warn = 1)
}


#' @rdname plotCovarianceMltpSingleLine
#' @name plotCovarianceComparisonMltpSingleLine
#' @export
#'
###############################################################################
#plotCovarianceComparisonMltpSingleLine
###############################################################################
plotCovarianceComparisonMltpSingleLine <- function(msa, helix1, helix2, comp.seq = FALSE, top.name = "FALSE",
                                                   sort = TRUE, arcs = TRUE, msa.all.seq = TRUE, conflict.cutoff = 0.01,
                                                   stable = FALSE, y = 0, flip = FALSE, grid = TRUE, dist.between = 0.1,
                                                   dist.part = 0.2,add = FALSE, pad = c(0, 0, 0, 0), shape = "circle",
                                                   debug = FALSE, col = "black",
                                                   base.colour = FALSE, palette = NA, text = TRUE, grid.col = "white",
                                                   grid.lwd = 0, text.cex = 0.2, text.col = "gray", text.font = 2,
                                                   text.family = "sans", species.cex = 0.2, species.col = "black",
                                                   species.font = 2, species.family = "mono", species = 10,
                                                   legend = TRUE,breaks = NA,cols = NA,log = FALSE,
                                                   include.lowest = TRUE,scale = FALSE,line = FALSE,arrow = FALSE,...){
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
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }



  dfTP <- merge(x = helix1, y = helix2, by = c("i", "j","length", "ent.i", "ent.j"))

  if ("value.x" %in% names(dfTP)) {
    if ("value.y" %in% names(dfTP)) {
      dfTP <- dfTP[, 1:6]
      colnames(dfTP) <- c("i", "j", "length",
                          "ent.i", "ent.j", "value")
    }
    else {
      colnames(dfTP) <- c("i", "j", "length",
                          "ent.i", "ent.j", "value")
    }
  }
  dfTP <- dfTP[, c("i", "j", "length", "value", "ent.i", "ent.j")]
  dfFN <- helix1[!dfTP, on = c("i", "j", "length", "ent.i", "ent.j")]
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


  #helix for covariance
  helix.cov <- rbind(helix1,helix2)
  attr(helix.cov, "length") <- attr(helix1, "length")

  #############################################################################################
  #update msa

  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }
  options(warn = -1)
  shape.conflict <- shape
  msa.data.table <- bstring.to.data.table(msa, comp.seq = comp.seq)
  uniq.fake <- copy(attr(helix.cov, "length"))
  if (length(unique(uniq.fake[, group])) != length(unique(msa.data.table[,group]))) {
    stop("Not equal amount of entities between fasta and msa")
  }
  if (!unique(sort(unique(msa.data.table[, group])) == sort(unique(msa.data.table[,group])))) {
    stop("Entities names are not same between fasta and msa")
  }
  uniq.fake <- uniq.fake[, `:=`(("width"), lapply(.SD, as.numeric)), .SDcols = "width"]
  index <- 1:uniq.fake[, .N]
  uniq.fake[, `:=`(N, index)]
  uniq.fake <- update.uniq.by.group.single(uniq.fake, top.name, sort)
  msa.group.seq.full <- list()
  seq.names.msa <- list()
  for (i in 1:nrow(uniq.fake)) {
    msa.group.seq.full[[uniq.fake[i, group]]] <- msa.data.table[group %like% uniq.fake[i, group]][, !"group"][, !"rest"]
    if (comp.seq) {
      seq.names.msa[[uniq.fake[i, group]]] <- msa.group.seq.full[[uniq.fake[i, group]]][, seq.name]
    }
    if (!comp.seq) {
      seq.names.msa[[uniq.fake[i, group]]] <- msa.group.seq.full[[uniq.fake[i,group]]][, name]
    }
  }
  for (i in 1:nrow(uniq.fake)) {
    if (length(unique(duplicated(seq.names.msa[[uniq.fake[i, group]]]))) > 1) {
      message("\nDuplicate sequence name present.\nJust first sequence will be used")
    }
    seq.names.msa[[uniq.fake[i, group]]] <- unique(seq.names.msa[[uniq.fake[i, group]]])
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
  }else {
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


  helix.trans.def <- copy(helix.cov)
  full.len.fake <- sum(uniq.fake[, width])
  for (i in 1:uniq.fake[, .N]) {
    uniq.fake[i, `:=`(start, full.len.fake - sum(uniq.fake$width[i:uniq.fake[, .N]]))]
  }
  uniq.fake[, `:=`(end, width + start)]
  helix.fake <- update.helix.coord.single(helix.trans.def,
                                          uniq.fake)
  helix.fake[, `:=`(i, i.coord)]
  helix.fake[, `:=`(j, j.coord)]
  helix.fake <- helix.fake[, `:=`(c("i", "j"), lapply(.SD, as.integer)), .SDcols = c("i", "j")]
  helix.fake[, `:=`(conflict, isConflictingHelix(helix.fake[, 1:4]))]
  helix.fake.restored <- copy(helix.fake)
  helix.fake.restored[, `:=`(i, helix.cov[, i])][, `:=`(j, helix.cov[, j])]
  helix.fake.restored <- helix.fake.restored[, !c("i.coord", "j.coord")]
  helix.conflict <- helix.fake.restored[conflict >= conflict.cutoff]
  helix.notconflict <- helix.fake.restored[conflict < conflict.cutoff]

  if(arcs){
    dfBP <- rbind(dfFN[, 1:6], dfTP[, 1:6])
    attr(dfBP, "length") <- attr(helix.cov, "length")
    max.height.helix.pos <- def.max.height.helix(dfBP, top.name = top.name,
                                                 sort = sort, dist.part = dist.part, stable = stable)
    max.height.helix.neg <- -def.max.height.helix(dfFP, top.name = top.name,
                                                  sort = sort, dist.part = dist.part, stable = stable)

    dist.between <- attr(max.height.helix.pos, "dist.between")
    uniq <- attr(max.height.helix.pos, "uniq")
    full.len <- uniq[.N, end]

  }else{
    max.height.helix.pos <- 1
    max.height.helix.pos <- -1
  }



  #define height of msa
  height.msa <- ifelse(grid, 1.25, 1)
  if (msa.all.seq) {
    msa.size <- height.msa * height.msa.max
  }else {
    msa.size <- height.msa * height.msa.min
  }


  if (!add) {
    if (species > 0) {
      pad[2] <- pad[2] + species
    }
    blankPlotMod(width = full.len, top = max.height.helix.pos+msa.size+1,
                 bottom = max.height.helix.neg-1, pad=pad, debug = debug,scale=scale,...)
  }


  if(length(unique(dfTP$value)) ==1 ){
    stop("True Positives have same value\n Colour to Blue")
    dfFP[,col := "#4169E1"]
  }else{
    dfFP <- colourByValueMltp(dfFP, get = TRUE, log = log,include.lowest = include.lowest)
  }
  if(length(unique(dfFP$value)) ==1 ){
    warning("False Positives have same value\n Colour to Blue")
    dfFP[,col := "#4169E1"]
  }else{
    dfTP <- colourByValueMltp(dfTP, get = TRUE, log = log,include.lowest = include.lowest)
  }



  if (arcs) {
    plotHelixMltpSingleLine(dfFN, col = "black", flip = FALSE, line = FALSE,
                            arrow = FALSE, add = TRUE, scale = scale, dist.part = dist.part,
                            stable = stable, shape = shape, debug = debug, top.name = top.name,
                            sort = sort,y=msa.size+1)
    plotHelixMltpSingleLine(dfTP, col = col, flip = FALSE, line = FALSE,
                            arrow = FALSE, add = TRUE, scale = scale, dist.part = dist.part,
                            stable = stable, shape = shape, debug = debug, top.name = top.name,
                            sort = sort,y=msa.size+1)
    plotHelixMltpSingleLine(dfFP, col = col, flip = TRUE, line = FALSE,
                            arrow = FALSE, add = TRUE, scale = scale, dist.part = dist.part,
                            stable = stable, shape = shape, debug = debug, top.name = top.name,
                            sort = sort,y = -1)
  }

  uniq$y <- 0
  if (line) {
    plot.bases.single(uniq, col.line)
  }
  if (arrow) {
    plot.arrows.single(uniq, dist.between, col.arrow)
  }


  options(warn = -1)
  uniq <- attr(helix.cov, "length")
  uniq <- update.uniq.by.group.single(uniq, top.name, sort)
  uniq <- prep.uniq.single(uniq, dist.part)

  if (base.colour) {
    cols <- getBaseColoursMod(msa.prep.equal)
  }else {
    cols <- getCovarianceColoursMltp(msa.prep.equal, helix.fake[,1:4])
  }
  if (all(is.na(palette))) {
    if (base.colour) {
		palette <- c("#E41A1C", "#4DAF4A", "#FF7F00", "#377EB8", "#BDBDBD","#984EA3","#00008b")
    }
    else {
      palette <- c("#00A651", "#0072BC", "#00B9F2",
                   "#F15A22", "#231F20", "#AAAAAA",
                   "#DA6FAB")
    }
  }else {
    if (length(palette) < ifelse(base.colour, 5, 7)) {
      palette <- rep(palette, length.out = ifelse(base.colour, 5, 7))
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
    if (is.vector(cols[, uniq.fake[i, start]:uniq.fake[i,
                                                       end]])) {
      msa.cols.storage[[i]] <- matrix(cols[, uniq.fake[i,
                                                       start]:uniq.fake[i, end]], nrow = 1)
    }
    else {
      msa.cols.storage[[i]] <- cols[, uniq.fake[i, start]:uniq.fake[i,
                                                                    end]]
    }
    if (nrow(msa.group.seq.full[[uniq.fake[i, group]]]) !=
        nrow(msa.cols.storage[[i]])) {
      for (j in 1:(nrow(msa.group.seq.full[[uniq.fake[i,
                                                      group]]]) - nrow(msa.cols.storage[[i]]))) {
        msa.cols.storage[[i]] <- rbind(msa.cols.storage[[i]],
                                       rep("#231F20", ncol(msa.cols.storage[[i]])))
      }
    }
    if (comp.seq) {
      tmp <- rbind(msa.selected.list[[uniq.fake[i, group]]],
                   msa.group.seq.full[[uniq.fake[i, group]]][!msa.selected.list[[uniq.fake[i,
                                                                                           group]]], on = c("seq.name", "sequence")])
    }
    else {
      tmp <- rbind(msa.selected.list[[uniq.fake[i, group]]],
                   msa.group.seq.full[[uniq.fake[i, group]]][!msa.selected.list[[uniq.fake[i,
                                                                                           group]]], on = c("name", "sequence")])
    }
    tmp1 <- strsplit(tmp[, sequence], "")
    matrix <- matrix(ncol = length(tmp1[[1]]), nrow = length(tmp1))
    for (j in 1:length(tmp1)) {
      matrix[j, ] <- tmp1[[j]]
    }
    msa.seq.storage[[i]] <- matrix
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



  for (i in 1:length(msa.cols.storage)) {
    for (j in 1:nrow(msa.seq.storage[[i]])) {
      if (grid) {
        plotCovarianceGrid(uniq[i, start] + 1, msa.size -height.msa * j, msa.cols.storage[[i]][j,], height = height.msa, border = grid.col,
                           lwd = grid.lwd)
        if (text) {
          plotCovarianceText(uniq[i, start] + 1, msa.size -height.msa * j, msa.seq.storage[[i]][j,], height = height.msa, cex = text.cex,
                             font = text.font, family = text.family,col = text.col)
        }
        if (species > 0) {
          plotCovarianceSpecies(species - uniq[i, start] +
                                  1, msa.size - height.msa * j, msa.name.storage[[i]][j],
                                height = height.msa, cex = species.cex,
                                font = species.font, family = species.family,
                                col = species.col)
        }
      }
      else {
        plotCovarianceLine(uniq[i, start] + 1, msa.size -j, msa.cols.storage[[i]][j, ])
      }
    }
  }


  if (legend) {
    if (base.colour) {
      text.legend <- c("A", "U", "G", "C", "-", "?","T")
    }
    else {
      text.legend <- c("Conservation", "Covariation",
                       "One-sided", "Invalid", "Unpaired",
                       "Gap", "Ambiguous")
    }
    legend("bottom", text.legend[used], border = NA,
           fill = palette[used], horiz = TRUE, bty = "n",
           xpd = NA)
  }
  options(warn = 1)

}


#' @rdname plotCovarianceMltpSingleLine
#' @name plotDoubleCovarianceComparisonMltpSingleLine
#' @export
#'
###############################################################################
#plotDoubleCovarianceComparisonMltpSingleLine
###############################################################################
plotDoubleCovarianceComparisonMltpSingleLine <- function(msa, helix1, helix2, comp.seq = FALSE, top.name = "FALSE",
                                                         sort = TRUE, arcs = TRUE, msa.all.seq = TRUE, conflict.cutoff = 0.01,
                                                         stable = FALSE, y = 0, flip = FALSE, grid = TRUE, dist.between = 0.1,
                                                         dist.part = 0.2, add = FALSE, pad = c(0, 0, 0, 0), shape = "circle",
                                                         debug = FALSE, col = "black",
                                                         base.colour = FALSE, palette = NA, text = TRUE, grid.col = "white",
                                                         grid.lwd = 0, text.cex = 0.2, text.col = "gray", text.font = 2,
                                                         text.family = "sans", species.cex = 0.2, species.col = "black",
                                                         species.font = 2, species.family = "mono", species = 10,
                                                         legend = TRUE, breaks = NA, cols = NA, log = FALSE, include.lowest = TRUE,
                                                         scale = FALSE, line = FALSE, arrow = FALSE,dist.y.between = 10,...){



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
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }

  options(warn = -1)

  #####################################################################

  # dfTP <- merge(x = helix1, y = helix2, by = c("i", "j", "length", "ent.i", "ent.j"))
  #
  # if ("value.x" %in% names(dfTP)) {
  #   if ("value.y" %in% names(dfTP)) {
  #     dfTP <- dfTP[, 1:6]
  #     colnames(dfTP) <- c("i", "j", "length", "ent.i", "ent.j", "value")
  #   }
  #   else {
  #     colnames(dfTP) <- c("i", "j", "length","ent.i", "ent.j", "value")
  #   }
  # }
  #
  # dfTP <- dfTP[, c("i", "j", "length", "value", "ent.i", "ent.j")]
  # dfFN <- helix1[!dfTP, on = c("i", "j", "length","ent.i", "ent.j")]
  # dfFP <- helix2[!dfTP, on = c("i", "j", "length", "ent.i", "ent.j")]
  #
  # if (nrow(dfTP) == 0) {
  #   stop("No True Positives")
  # }
  # if (nrow(dfFN) == 0) {
  #   stop("No False Negatives")
  # }
  # if (nrow(dfFP) == 0) {
  #   stop("No False Positives")
  # }
  # attr(dfTP, "length") <- attr(helix1, "length")
  # attr(dfFN, "length") <- attr(helix1, "length")
  # attr(dfFN, "length") <- attr(helix1, "length")
  #
  # if(unique(dfTP$value) ==1 ){
  #   warning("True Positives have same value\n Colour to Blue")
  #   dfTP[,col := "#4169E1"]
  # }else{
  #   if(is.na(breaks)){
  #     if(is.na(cols)){
  #       dfTP <- colourByValueMltp(dfTP,get = TRUE, log = log, include.lowest = include.lowest)
  #     }else{
  #       dfTP <- colourByValueMltp(dfTP,get = TRUE, cols = cols, log = log, include.lowest = include.lowest)
  #     }
  #   }else{
  #     if(is.na(cols)){
  #       dfTP <- colourByValueMltp(dfTP,get = TRUE,breaks = breaks, log = log, include.lowest = include.lowest)
  #     }else{
  #       dfTP <- colourByValueMltp(dfTP,get = TRUE,breaks = breaks, cols = cols, log = log, include.lowest = include.lowest)
  #     }
  #   }
  #
  # }
  # if(unique(dfFP$value) ==1 ){
  #   warning("False Positives have same value\n Colour to Blue")
  #   dfFP[,col := "#4169E1"]
  # }else{
  #   if(is.na(breaks)){
  #     if(is.na(cols)){
  #       dfFP <- colourByValueMltp(dfFP,get = TRUE, log = log, include.lowest = include.lowest)
  #     }else{
  #       dfFP <- colourByValueMltp(dfFP,get = TRUE, cols = cols, log = log, include.lowest = include.lowest)
  #     }
  #   }else{
  #     if(is.na(cols)){
  #       dfFP <- colourByValueMltp(dfFP,get = TRUE,breaks = breaks, log = log, include.lowest = include.lowest)
  #     }else{
  #       dfFP <- colourByValueMltp(dfFP,get = TRUE,breaks = breaks, cols = cols, log = log, include.lowest = include.lowest)
  #     }
  #   }
  # }
  #
  # #left side dfTP(coloured) and dfFN(full black)
  # dfFN$col <- "#000000"
  # df.top <- rbind(dfTP,dfFN)
  # attr(df.top, "length") <- attr(helix1, "length")
  #
  # #right side dfFP
  # df.bottom <- dfFP
  # attr(df.bottom, "length") <- attr(helix1, "length")
  #
  #
  # #################################################################################
  #
  # if (arcs) {
  #   max.height.helix.pos <- def.max.height.helix(df.top, top.name = top.name,sort = sort, dist.part = dist.part, stable = stable)
  #   max.height.helix.neg <- -def.max.height.helix(df.bottom, top.name = top.name, sort = sort, dist.part = dist.part, stable = stable)
  #   dist.between <- attr(max.height.helix.pos, "dist.between")
  #   uniq <- attr(max.height.helix.pos, "uniq")
  #   full.len <- uniq[.N, end]
  # }else {
  #   max.height.helix.pos <- 1
  #   max.height.helix.neg <- -1
  # }
  #
  # #################################################################################
  #
  # #####Extract sequences and make msa as 1 BStringSet
  # #####For next using on default functions
  # #transfor fasta to data.table format
  # msa.data.table <- bstring.to.data.table(msa,comp.seq = comp.seq)
  # #call uniq from helix file
  # uniq.fake <- copy(attr(helix1,"length"))
  #
  # #check if N of unique entities are same
  # if(length(unique(uniq.fake[,group])) != length(unique(msa.data.table[,group]))){
  #   stop("Not equal amount of entities between fasta and msa")
  # }
  # if(!unique(sort(unique(msa.data.table[,group])) == sort(unique(msa.data.table[,group])))){
  #   stop("Entities names are not same between fasta and msa")
  # }
  #
  # uniq.fake <- uniq.fake[,("width"):= lapply(.SD, as.numeric), .SDcols = "width"]
  # index <- 1:uniq.fake[,.N]
  # uniq.fake[,N:=index]
  # #update uniq if it needed by top.name
  # uniq.fake <- update.uniq.by.group.single(uniq.fake,top.name,sort)
  #
  # #contain full msa sequences
  # msa.group.seq.full <- list()
  # #contain names of sequences
  # seq.names.msa <- list()
  # for(i in 1:nrow(uniq.fake)){
  #   msa.group.seq.full[[uniq.fake[i,group]]] <- msa.data.table[group %like% uniq.fake[i,group]][,!"group"][,!"rest"]
  #   if(comp.seq){
  #     seq.names.msa[[uniq.fake[i,group]]] <- msa.group.seq.full[[uniq.fake[i,group]]][,seq.name]
  #   }
  #   if(!comp.seq){
  #     seq.names.msa[[uniq.fake[i,group]]] <- msa.group.seq.full[[uniq.fake[i,group]]][,name]
  #   }
  #
  # }
  # #define duplicates and take massage to user, after delete duplicate names
  # for(i in 1:nrow(uniq.fake)){
  #   if(length(unique(duplicated(seq.names.msa[[uniq.fake[i,group]]])))>1){
  #     message("\nDuplicate sequence name present.\nJust first sequence will be used")
  #   }
  #   seq.names.msa[[uniq.fake[i,group]]] <- unique(seq.names.msa[[uniq.fake[i,group]]])
  # }
  #
  # if(comp.seq){
  #   #find common names between all entities
  #   common.seq.names <- Reduce(intersect, seq.names.msa)
  #   if(length(common.seq.names)==0){
  #     warning("No similar seq names present in entities\nFirst will be choosen")
  #   }
  # }
  # #define minimum of msa seq
  # msa.rows.height <- c()
  # for(i in 1:nrow(uniq.fake)){
  #   msa.rows.height <-c(msa.rows.height,nrow(msa.group.seq.full[[uniq.fake[i,group]]]))
  # }
  #
  # height.msa.min <- min(msa.rows.height)
  # height.msa.max <- max(msa.rows.height)
  #
  # height.msa <- ifelse(grid, 1.25, 1)
  # if (msa.all.seq) {
  #   msa.size <- height.msa * height.msa.max
  # }else {
  #   msa.size <- height.msa * height.msa.min
  # }
  #
  #
  # ######################################################
  # if (!add) {
  #   if (species > 0) {
  #     pad[2] <- pad[2] + species
  #   }
  #   blankPlotMod(width = full.len, top = max.height.helix.pos +msa.size + dist.y.between + 1,
  #                bottom = max.height.helix.neg - msa.size - dist.y.between - 1,
  #                pad = pad, debug = debug,scale=scale,...)
  # }
  #
  #
  # plotCovarianceMltpSingleLine(msa = msa, helix = df.top, comp.seq = comp.seq, top.name = top.name,
  #                              sort = sort, arcs = arcs, msa.all.seq = msa.all.seq, conflict.cutoff = conflict.cutoff,
  #                              stable = stable, y = dist.y.between + 1, flip = FALSE, grid = grid, dist.between = dist.between,
  #                              dist.part = dist.part, add = TRUE, pad = pad, shape = shape,
  #                              debug = debug, col = col, conflict.col = conflict.col, conflict.lty = conflict.lty,
  #                              conflict.lwd = conflict.lwd, base.colour = base.colour, palette = palette, text = text,
  #                              grid.col = grid.col, grid.lwd = grid.lwd, text.cex = text.cex, text.col = text.col,
  #                              text.font = text.font, text.family = text.family, species.cex = species.cex,
  #                              species.col = species.col, species.font = species.font, species.family = species.family,
  #                              species = species, legend = FALSE)
  #
  # plotCovarianceMltpSingleLine(msa = msa, helix = df.bottom, comp.seq = comp.seq, top.name = top.name,
  #                              sort = sort, arcs = arcs, msa.all.seq = msa.all.seq, conflict.cutoff = conflict.cutoff,
  #                              stable = stable, y = -dist.y.between - 1, flip = TRUE, grid = grid, dist.between = dist.between,
  #                              dist.part = dist.part, add = TRUE, pad = pad, shape = shape,
  #                              debug = debug, col = col, conflict.col = conflict.col, conflict.lty = conflict.lty,
  #                              conflict.lwd = conflict.lwd, base.colour = base.colour, palette = palette, text = text,
  #                              grid.col = grid.col, grid.lwd = grid.lwd, text.cex = text.cex, text.col = text.col,
  #                              text.font = text.font, text.family = text.family, species.cex = species.cex,
  #                              species.col = species.col, species.font = species.font, species.family = species.family,
  #                              species = species, legend = FALSE)
  #
  #
  # if (all(is.na(palette))) {
  #   if (base.colour) {
  #     palette <- c("#E41A1C", "#4DAF4A", "#FF7F00", "#377EB8", "#BDBDBD","#984EA3","#00008b")
  #   }else {
  #     palette <- c("#00A651", "#0072BC", "#00B9F2",
  #                  "#F15A22", "#231F20", "#AAAAAA",
  #                  "#DA6FAB")
  #   }
  # }else {
  #   if (length(palette) < ifelse(base.colour, 5, 7)) {
  #     palette <- rep(palette, length.out = ifelse(base.colour,5, 7))
  #   }
  # }
  # if (legend) {
  #   if (base.colour) {
  #     text.legend <- c("A", "U", "G","C", "-", "?","T")
  #   }else {
  #     text.legend <- c("Conservation", "Covariation",
  #                      "One-sided", "Invalid", "Unpaired",
  #                      "Gap", "Ambiguous")
  #   }
  #   legend("bottom", text.legend, border = NA,
  #          fill = palette, horiz = TRUE, bty = "n", xpd = NA)
  # }
  #
  # options(warn = 1)

  #####################################################################

  dfTP <- merge(x = helix1, y = helix2, by = c("i", "j","length", "ent.i", "ent.j"))

  if ("value.x" %in% names(dfTP)) {
    if ("value.y" %in% names(dfTP)) {
      dfTP <- dfTP[, 1:6]
      colnames(dfTP) <- c("i", "j", "length",
                          "ent.i", "ent.j", "value")
    }else {
      colnames(dfTP) <- c("i", "j", "length",
                          "ent.i", "ent.j", "value")
    }
  }
  dfTP <- dfTP[, c("i", "j", "length", "value",
                   "ent.i", "ent.j")]
  dfFN <- helix1[!dfTP, on = c("i", "j", "length",
                               "ent.i", "ent.j")]
  dfFP <- helix2[!dfTP, on = c("i", "j", "length",
                               "ent.i", "ent.j")]
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

  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }

  #helices for top and down
  helix.cov1 <- copy(helix1)
  attr(helix.cov1, "length") <- attr(helix1, "length")
  helix.cov2 <- copy(helix2)
  attr(helix.cov2, "length") <- attr(helix2, "length")


  #create msa for covariance
  shape.conflict <- shape
  msa.data.table <- bstring.to.data.table(msa, comp.seq = comp.seq)
  uniq.fake <- copy(attr(helix.cov1, "length"))
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
    seq.names.msa[[uniq.fake[i, group]]] <- unique(seq.names.msa[[uniq.fake[i,
                                                                            group]]])
  }
  if (comp.seq) {
    common.seq.names <- Reduce(intersect, seq.names.msa)
    if (length(common.seq.names) == 0) {
      warning("No similar seq names present in entities\nFirst will be choosen")
    }
  }
  msa.rows.height <- c()
  for (i in 1:nrow(uniq.fake)) {
    msa.rows.height <- c(msa.rows.height, nrow(msa.group.seq.full[[uniq.fake[i,
                                                                             group]]]))
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
  }else {
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


  #make colours matrices for msa

  helix.trans.def <- copy(helix.cov1)
  uniq1.fake <- copy(attr(helix.cov1, "length"))

  full.len.fake <- sum(uniq1.fake$width)

  for (i in 1:uniq1.fake[, .N]) {
    uniq1.fake[i, `:=`(start, full.len.fake - sum(uniq.fake$width[i:uniq1.fake[,.N]]))]
  }
  uniq1.fake[, `:=`(end, width + start)]
  helix1.fake <- update.helix.coord.single(helix.trans.def,uniq1.fake)
  helix1.fake[, `:=`(i, i.coord)]
  helix1.fake[, `:=`(j, j.coord)]
  helix1.fake <- helix1.fake[, `:=`(c("i", "j"), lapply(.SD, as.integer)), .SDcols = c("i", "j")]


  helix.trans.def <- copy(helix.cov2)
  uniq2.fake <- copy(attr(helix.cov2, "length"))

  full.len.fake <- sum(uniq2.fake[, width])
  for (i in 1:uniq2.fake[, .N]) {
    uniq2.fake[i, `:=`(start, full.len.fake - sum(uniq.fake$width[i:uniq2.fake[,.N]]))]
  }
  uniq2.fake[, `:=`(end, width + start)]
  helix2.fake <- update.helix.coord.single(helix.trans.def,uniq2.fake)
  helix2.fake[, `:=`(i, i.coord)]
  helix2.fake[, `:=`(j, j.coord)]
  helix2.fake <- helix2.fake[, `:=`(c("i", "j"), lapply(.SD, as.integer)), .SDcols = c("i", "j")]



  if (base.colour) {
    cols1 <- getBaseColoursMod(msa.prep.equal)
    cols2 <- getBaseColoursMod(msa.prep.equal)
  }else {
    cols1 <- getCovarianceColoursMltp(msa.prep.equal, helix1.fake[,1:4])
    cols2 <- getCovarianceColoursMltp(msa.prep.equal, helix2.fake[,1:4])
  }

  #assign colours
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



  ############################################################################
  ############################################################################
  ############################################################################
  used1 <- sort(unique(c(cols1)))
  dim1 <- dim(cols1)
  cols1 <- palette[cols1]
  dim(cols1) <- dim1
  uniq1.fake[, `:=`(start, start + 1)]
  msa.cols.storage.1 <- list()
  msa.name.storage.1 <- list()
  msa.seq.storage.1 <- list()


  for (i in 1:uniq1.fake[, .N]) {
    if (is.vector(cols1[, uniq1.fake[i, start]:uniq1.fake[i,end]])) {
      msa.cols.storage.1[[i]] <- matrix(cols1[, uniq1.fake[i,start]:uniq1.fake[i, end]], nrow = 1)
    }else {
      msa.cols.storage.1[[i]] <- cols1[, uniq1.fake[i, start]:uniq1.fake[i, end]]
    }
    if (nrow(msa.group.seq.full[[uniq1.fake[i, group]]]) != nrow(msa.cols.storage.1[[i]])) {
      for (j in 1:(nrow(msa.group.seq.full[[uniq1.fake[i, group]]]) - nrow(msa.cols.storage.1[[i]]))) {
        msa.cols.storage.1[[i]] <- rbind(msa.cols.storage.1[[i]], rep("#231F20", ncol(msa.cols.storage.1[[i]])))
      }
    }
    if (comp.seq) {
      tmp <- rbind(msa.selected.list[[uniq1.fake[i, group]]],msa.group.seq.full[[uniq1.fake[i, group]]]
                   [!msa.selected.list[[uniq1.fake[i,group]]],on = c("seq.name", "sequence")])
    }else {
      tmp <- rbind(msa.selected.list[[uniq1.fake[i, group]]], msa.group.seq.full[[uniq1.fake[i, group]]]
                   [!msa.selected.list[[uniq1.fake[i,group]]], on = c("name", "sequence")])
    }
    tmp1 <- strsplit(tmp[, sequence], "")
    matrix <- matrix(ncol = length(tmp1[[1]]), nrow = length(tmp1))
    for (j in 1:length(tmp1)) {
      matrix[j, ] <- tmp1[[j]]
    }
    msa.seq.storage.1[[i]] <- matrix
    if (comp.seq) {
      msa.name.storage.1[[i]] <- tmp[, seq.name]
    }else {
      msa.name.storage.1[[i]] <- tmp[, name]
    }
  }


  if (!msa.all.seq) {
    for (i in 1:length(msa.cols.storage.1)) {
      msa.seq.storage.1[[i]] <- msa.seq.storage.1[[i]][1:height.msa.min,]
      msa.cols.storage.1[[i]] <- msa.cols.storage.1[[i]][1:height.msa.min,]
      msa.name.storage.1[[i]] <- msa.name.storage.1[[i]][1:height.msa.min]
    }
  }


  ############################################################################
  ############################################################################
  ############################################################################

  used2 <- sort(unique(c(cols2)))
  dim2 <- dim(cols2)
  cols2 <- palette[cols2]
  dim(cols2) <- dim2
  uniq2.fake[, `:=`(start, start + 1)]
  msa.cols.storage.2 <- list()
  msa.name.storage.2 <- list()
  msa.seq.storage.2 <- list()


  for (i in 1:uniq2.fake[, .N]) {
    if (is.vector(cols2[, uniq2.fake[i, start]:uniq2.fake[i,end]])) {
      msa.cols.storage.2[[i]] <- matrix(cols2[, uniq2.fake[i,start]:uniq2.fake[i, end]], nrow = 1)
    }else {
      msa.cols.storage.2[[i]] <- cols2[, uniq2.fake[i, start]:uniq2.fake[i, end]]
    }
    if (nrow(msa.group.seq.full[[uniq2.fake[i, group]]]) != nrow(msa.cols.storage.2[[i]])) {
      for (j in 1:(nrow(msa.group.seq.full[[uniq2.fake[i, group]]]) - nrow(msa.cols.storage.2[[i]]))) {
        msa.cols.storage.2[[i]] <- rbind(msa.cols.storage.2[[i]], rep("#231F20", ncol(msa.cols.storage.2[[i]])))
      }
    }
    if (comp.seq) {
      tmp <- rbind(msa.selected.list[[uniq2.fake[i, group]]],msa.group.seq.full[[uniq2.fake[i, group]]]
                   [!msa.selected.list[[uniq2.fake[i,group]]],on = c("seq.name", "sequence")])
    }else {
      tmp <- rbind(msa.selected.list[[uniq2.fake[i, group]]], msa.group.seq.full[[uniq2.fake[i, group]]]
                   [!msa.selected.list[[uniq2.fake[i,group]]], on = c("name", "sequence")])
    }
    tmp1 <- strsplit(tmp[, sequence], "")
    matrix <- matrix(ncol = length(tmp1[[1]]), nrow = length(tmp1))
    for (j in 1:length(tmp1)) {
      matrix[j, ] <- tmp1[[j]]
    }
    msa.seq.storage.2[[i]] <- matrix
    if (comp.seq) {
      msa.name.storage.2[[i]] <- tmp[, seq.name]
    }else {
      msa.name.storage.2[[i]] <- tmp[, name]
    }
  }


  if (!msa.all.seq) {
    for (i in 1:length(msa.cols.storage.2)) {
      msa.seq.storage.2[[i]] <- msa.seq.storage.2[[i]][1:height.msa.min,]
      msa.cols.storage.2[[i]] <- msa.cols.storage.2[[i]][1:height.msa.min,]
      msa.name.storage.2[[i]] <- msa.name.storage.2[[i]][1:height.msa.min]
    }
  }


  #################################################################################

  if (arcs) {
    dfBP <- rbind(dfFN[, 1:6], dfTP[, 1:6])
    attr(dfBP, "length") <- attr(helix.cov1, "length")
    max.height.helix.pos <- def.max.height.helix(dfBP, top.name = top.name,sort = sort, dist.part = dist.part, stable = stable)
    max.height.helix.neg <- -def.max.height.helix(dfFP, top.name = top.name, sort = sort, dist.part = dist.part, stable = stable)
    dist.between <- attr(max.height.helix.pos, "dist.between")
    uniq <- attr(max.height.helix.pos, "uniq")
    full.len <- uniq[.N, end]
  }else {
    max.height.helix.pos <- 1
    max.height.helix.neg <- -1
  }

  height.msa <- ifelse(grid, 1.25, 1)
  if (msa.all.seq) {
    msa.size <- height.msa * height.msa.max
  }else {
    msa.size <- height.msa * height.msa.min
  }

  ######################################################
  if (!add) {
    if (species > 0) {
      pad[2] <- pad[2] + species
    }
    blankPlotMod(width = full.len, top = max.height.helix.pos +
                   msa.size + dist.y.between + 1, bottom = max.height.helix.neg - msa.size - dist.y.between - 1,
                 pad = pad, debug = debug,scale=scale,...)
  }

  if(length(unique(dfTP$value)) ==1 ){
    stop("True Positives have same value\n Colour to Blue")
    dfFP[,col := "#4169E1"]
  }else{
    dfFP <- colourByValueMltp(dfFP, get = TRUE, log = log,include.lowest = include.lowest)
  }
  if(length(unique(dfFP$value)) ==1 ){
    warning("False Positives have same value\n Colour to Blue")
    dfTP[,col := "#4169E1"]
  }else{
    dfTP <- colourByValueMltp(dfTP, get = TRUE, log = log,include.lowest = include.lowest)
  }

  if (arcs) {
    plotHelixMltpSingleLine(dfFN, col = "black", flip = FALSE,
                            line = FALSE, arrow = FALSE, add = TRUE, scale = scale,
                            dist.part = dist.part, stable = stable, shape = shape,
                            debug = debug, top.name = top.name, sort = sort,
                            y = dist.y.between + msa.size + 1)
    plotHelixMltpSingleLine(dfTP, col = col, flip = FALSE,
                            line = FALSE, arrow = FALSE, add = TRUE, scale = scale,
                            dist.part = dist.part, stable = stable, shape = shape,
                            debug = debug, top.name = top.name, sort = sort,
                            y = dist.y.between + msa.size + 1)
    plotHelixMltpSingleLine(dfFP, col = col, flip = TRUE,
                            line = FALSE, arrow = FALSE, add = TRUE, scale = scale,
                            dist.part = dist.part, stable = stable, shape = shape,
                            debug = debug, top.name = top.name, sort = sort,
                            y = -dist.y.between-msa.size-1)
  }

  uniq$y <- 0
  if (line) {
    plot.bases.single(uniq, col.line)
  }
  if (arrow) {
    plot.arrows.single(uniq, dist.between, col.arrow)
  }



  for (i in 1:length(msa.cols.storage.1)) {
    for (j in 1:nrow(msa.seq.storage.1[[i]])) {
      if (grid) {
        plotCovarianceGrid(uniq[i, start] + 1, dist.y.between + msa.size - height.msa * j, msa.cols.storage.1[[i]][j, ], height = height.msa, border = grid.col, lwd = grid.lwd)
        plotCovarianceGrid(uniq[i, start] + 1, - dist.y.between - msa.size + height.msa * j, msa.cols.storage.2[[i]][j, ], height = height.msa, border = grid.col, lwd = grid.lwd)
        if (text) {
          plotCovarianceText(uniq[i, start] + 1, dist.y.between +msa.size -
                               height.msa * j, msa.seq.storage.1[[i]][j, ],
                             height = height.msa, cex = text.cex, font = text.font,
                             family = text.family, col = text.col)
          plotCovarianceText(uniq[i, start] + 1,- dist.y.between - msa.size +
                               height.msa * j, msa.seq.storage.2[[i]][j, ],
                             height = height.msa, cex = text.cex, font = text.font,
                             family = text.family, col = text.col)

        }
        if (species > 0) {
          plotCovarianceSpecies(species - uniq[i, start] + 1, dist.y.between + msa.size - height.msa * j, msa.name.storage.1[[i]][j],
                                height = height.msa, cex = species.cex, font = species.font,
                                family = species.family, col = species.col)
          plotCovarianceSpecies(species - uniq[i, start] + 1, - dist.y.between - msa.size + height.msa * j, msa.name.storage.2[[i]][j],
                                height = height.msa, cex = species.cex, font = species.font,
                                family = species.family, col = species.col)

        }
      }else {
        plotCovarianceLine(uniq[i, start] + 1, dist.y.between + msa.size - j, msa.cols.storage.1[[i]][j, ])
        plotCovarianceLine(uniq[i, start] + 1, - dist.y.between - msa.size + j, msa.cols.storage.2[[i]][j, ])
      }
    }
  }

  if (legend) {
    if (base.colour) {
      text.legend <- c("A", "U", "G", "C", "-", "?","T")
    }else {
      text.legend <- c("Conservation", "Covariation",
                       "One-sided", "Invalid", "Unpaired",
                       "Gap", "Ambiguous")
    }
    legend("bottom", text.legend[unique(c(used1,used2))], border = NA,
           fill = palette[unique(c(used1,used2))], horiz = TRUE, bty = "n",
           xpd = NA)
  }

  options(warn = 1)

}



#' @rdname plotCovarianceMltpSingleLine
#' @name plotCovarianceDoubleMltpSingleLine
#' @export
#'
###############################################################################
#plotDoubleHelixMltpSingleLine
###############################################################################
plotCovarianceDoubleMltpSingleLine <- function(helix1,helix2,msa1,msa2, comp.seq = FALSE, top.name = "FALSE",
                                               sort = FALSE, arcs = TRUE, msa.all.seq = TRUE, conflict.cutoff = 0.01,
                                               stable = FALSE, y = 0, flip = FALSE, grid = TRUE, dist.between = 0.1,
                                               dist.part = 0.2, add = FALSE, pad = c(0, 0, 0, 0), shape = "circle",
                                               debug = FALSE, col = "black", conflict.col = NA, conflict.lty = 2,conflict.lwd = 1,
                                               base.colour = FALSE, palette = NA, text = TRUE, grid.col = "white",
                                               grid.lwd = 0, text.cex = 0.2, text.col = "gray", text.font = 2,
                                               text.family = "sans", species.cex = 0.2, species.col = "black",
                                               species.font = 2, species.family = "mono", species = 10,
                                               legend = TRUE,dist.plot = 10,...){

  if (!is.helix.mltp(helix1) | !is.helix.mltp(helix2)) {
    stop("invalid helix input")
  }
  if (is(msa1, "XStringSet")) {
    msa.pre.1 <- as.character(msa1)
  }else {
    stop("invalid msa sturucture")
  }
  if (is(msa2, "XStringSet")) {
    msa.pre.2 <- as.character(msa2)
  }else {
    stop("invalid msa sturucture")
  }
  msa.data.table.1 <- bstring.to.data.table(msa.pre.1, comp.seq = comp.seq)
  msa.data.table.2 <- bstring.to.data.table(msa.pre.2, comp.seq = comp.seq)
  uniq1.fake <- copy(attr(helix1, "length"))
  uniq1.fake <- uniq1.fake[, `:=`(("width"), lapply(.SD,
                                                    as.numeric)), .SDcols = "width"]
  index <- 1:uniq1.fake[, .N]
  uniq1.fake[, `:=`(N, index)]
  uniq1.fake <- update.uniq.by.group.single(uniq1.fake, top.name,
                                            sort)
  uniq2.fake <- copy(attr(helix2, "length"))
  uniq2.fake <- uniq2.fake[, `:=`(("width"), lapply(.SD,
                                                    as.numeric)), .SDcols = "width"]
  index <- 1:uniq2.fake[, .N]
  uniq2.fake[, `:=`(N, index)]
  uniq2.fake <- update.uniq.by.group.single(uniq2.fake, top.name, sort)


  uniq1.fake <- prep.uniq.single(uniq1.fake,dist.part)
  uniq2.fake <- prep.uniq.single(uniq2.fake,dist.part)


  width.plot <- sum(min(uniq1.fake$start), max(uniq1.fake$end))
  msa1.uniq.group <- unique(msa.data.table.1$group)
  msa2.uniq.group <- unique(msa.data.table.2$group)
  msa1.ent.nrow <- c()
  msa2.ent.nrow <- c()
  for (k in 1:length(msa1.uniq.group)) {
    msa1.ent.nrow <- c(msa1.ent.nrow, nrow(msa.data.table.1[group == msa1.uniq.group[k]]))
  }
  for (k in 1:length(msa2.uniq.group)) {
    msa2.ent.nrow <- c(msa2.ent.nrow, nrow(msa.data.table.1[group == msa2.uniq.group[k]]))
  }
  height.msa <- ifelse(grid, 1.25, 1)
  if (msa.all.seq) {
    msa.pos <- max(msa1.ent.nrow) * height.msa
    msa.neg <- -max(msa2.ent.nrow) * height.msa
  }else {
    msa.pos <- min(msa1.ent.nrow) * height.msa
    msa.neg <- -min(msa2.ent.nrow) * height.msa
  }
  helix.pos <- def.max.height.helix(helix1, top.name = top.name,
                                    sort = sort, dist.part = dist.part, stable = stable)
  helix.neg <- -def.max.height.helix(helix2, top.name = top.name,
                                     sort = sort, dist.part = dist.part, stable = stable)
  if (!add) {
    if (species > 0) {
      pad[2] <- pad[2] + species
    }
    blankPlotMod(width = width.plot, top = helix.pos + msa.pos +
                   dist.plot/2, bottom = helix.neg + msa.neg - dist.plot/2,
                 pad, debug = debug, ...)
  }


  plotCovarianceMltpSingleLine(helix = helix1, msa = msa1,
                               flip = FALSE, add = TRUE, y = msa.pos + dist.plot/2,
                               comp.seq = comp.seq, top.name = top.name, sort = sort,
                               arcs = arcs, msa.all.seq = msa.all.seq, conflict.cutoff = conflict.cutoff,
                               stable = stable, grid = grid, dist.between = dist.between,
                               dist.part = dist.part, pad = pad, shape = shape, debug = debug,
                               col = col, conflict.col = conflict.col, conflict.lty = conflict.lty,conflict.lwd = conflict.lwd,
                               base.colour = base.colour, palette = palette, text = text,
                               grid.col = grid.col, grid.lwd = grid.lwd, text.cex = text.cex,
                               text.col = text.col, text.font = text.font, text.family = text.family,
                               species.cex = species.cex, species.col = species.col,
                               species.font = species.font, species.family = species.family,
                               species = species, legend = FALSE)
  plotCovarianceMltpSingleLine(helix = helix2, msa = msa2,
                               flip = TRUE, add = TRUE, y = -dist.plot/2, comp.seq = comp.seq,
                               top.name = top.name, sort = sort, arcs = arcs, msa.all.seq = msa.all.seq,
                               conflict.cutoff = conflict.cutoff, stable = stable, grid = grid,
                               dist.between = dist.between, dist.part = dist.part, pad = pad,
                               shape = shape, debug = debug, col = col, conflict.col = conflict.col,
                               conflict.lty = conflict.lty,conflict.lwd = conflict.lwd, base.colour = base.colour,
                               palette = palette, text = text, grid.col = grid.col,
                               grid.lwd = grid.lwd, text.cex = text.cex, text.col = text.col,
                               text.font = text.font, text.family = text.family, species.cex = species.cex,
                               species.col = species.col, species.font = species.font,
                               species.family = species.family, species = species, legend = legend)


}











