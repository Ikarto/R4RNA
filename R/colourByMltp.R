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




#' Assign colours to helices
#'
#' @description
#' Functions to generate colours for helices by various rules,
#' including integer counts, value ranges, percent identity covariation,
#' conservation, percentage canonical basepair, basepair frequency,
#' and non-pseudoknotted groups.
#'
#'
#' @param helix A helix data table to be coloured.
#' @param cols An array of characters (or numbers) representing a set of colours to
#' colour helix with. When missing, a default set of colours from defaultPalette()
#' will be used. Valid input include hex codes, colour names from the colours function,
#' and integer numbers. The colours will be interpreted as being from best to worst.
#' @param counts An array of integers the same length as cols, dictating the number
#' of times each corresponding colour should be used. When missing, the function will
#' divide the number of helices evenly over each of the colours available.
#' @param breaks An integer number of intervals to break the ‘value’ column of helix
#' into, or a list of numbers defining the interval breaks. If missing, the range of
#' ‘helix$value’ will automatically be split evenly into intervals for each colour available.
#' @param get If TRUE, returns the input helix with a col column, else simply
#' returns an array of colours the same length as the number of row in helix.
#' The exceptions are colourByBasepairFrequency and colourByUnknottedGroupsMltp
#' which will return a different helix if TRUE, and a list of colours that
#' will not match the input helix if FALSE.
#' @param log If TRUE, will breaks values into even log10 space intervals, useful when values are p-values.
#' @param include.lowest Whether the lowest interval should include the lowest value, passed to cut
#' @param top.name Name of entity, which will be set as top name and will refer to first entity in plot
#' @param sort Sort entities by alphabetic order
#' @param msa A multiple sequence alignment, such as those returned by readBstringSet
#' @param comp.seq Option to compare sequences from different entities by their "name" and provide
#' proper covariance colouring in future
#'
#' @export
#'
###############################################################################
#Assign color based on p-value
###############################################################################
colourByValueMltp <- function(helix,cols,
                              breaks,get = FALSE,
                              log = FALSE,include.lowest = TRUE){

  options(warn = -1)

  if (!is.helix.mltp(helix)) {
    stop("Not a valid helix data frame, aborting")
  }else{
    if (nrow(helix) == 0) {
      warning("No helices detected")
      if (get) {
        return(helix)
      }
      else {
        return(NULL)
      }
    }
    if (all(is.na(helix$value))) {
      stop("Cowardly refusing to deal with NA/NaN values")
    }
  }
  if (missing(cols)) {
    cols <- defaultPalette()
  }
  if (missing(breaks)) {
    breaks <- seq(min(helix[,value], na.rm = TRUE), max(helix[,value], na.rm = TRUE),
                  length.out = length(cols) + 1)
  }else {
    if (length(breaks) == 1 && breaks == 1) {
      breaks <- c(min(helix[,value], na.rm = TRUE), max(helix[,value],na.rm = TRUE))
    }
  }
  if (log) {
    test <- helix[value>0,value]
    start <- logfloor(min(test, na.rm = TRUE))
    end <- logceiling(max(test, na.rm = TRUE))
    breaks <- c(min(helix[,value], na.rm = TRUE), logseq(start, end)[-1])
  }

  levels <- cut(helix[,value], breaks = breaks, labels = NULL,
                ordered_result = FALSE, include.lowest = include.lowest)

  legend <- levels(levels)
  if (length(legend) > length(cols)) {
    warning(paste("Number of intervals (", length(legend),
                  ") greater than colours (", length(cols), "), some intervals will be colourless",
                  sep = ""))
  }
  fill <- cols[1:length(legend)]
  output <- cols[as.numeric(levels)]
  attr(output, "legend") <- legend
  attr(output, "fill") <- fill
  if (get) {
    helix[,col := output]
    attr(helix, "legend") <- legend
    attr(helix, "fill") <- fill
    return(helix)
  }else {
    return(output)
  }

  options(warn = 1)
}



###############################################################################
#unknottedGroupsMltp
###############################################################################
unknottedGroupsMltp <- function(helix) {
  if (!is.helix.mltp(helix)) {
    stop("Invalid input")
  }
  if (any(helix$length > 1)) {
    warning("Expanding helix to basepairs")
    helix <- copy(expandHelixMltp(helix))
  }
  if (nrow(helix) == 0) {
    stop("No basepairs in helix")
  }
  group <- rep(0, nrow(helix))
  group[1] <- 1
  if (nrow(helix) == 1) {
    return(group)
  }
  for (i in 2:nrow(helix)) {
    test <- 1
    while (any(is_pseudoknotted(helix[i, ][,1:4], helix[which(group == test), ][,1:4]))) {
      test <- test + 1
    }
    group[i] <- test
  }
  return(group)
}





#' @rdname colourByValueMltp
#' @export
#'
###############################################################################
#colourByUnknottedGroupsMltp
###############################################################################
colourByUnknottedGroupsMltp <- function(helix, cols, get = FALSE){
  if (!is.helix.mltp(helix)) {
    stop("Not a valid helix data frame, aborting")
  }
  if (any(helix$length >1)) {
    stop("Helix is not expanded, provide expand version")
  }

  expanded <- expandHelixMltp(helix)
  group <- as.factor(unknottedGroupsMltp(expanded))
  if (missing(cols)) {
    cols <- hcl.colors(max(as.integer(group)),palette = "dark2")
  }
  if (length(levels(group)) > length(cols)) {
    warning(paste("Number of groups (", length(levels(group)),
                  ") greater than colours (", length(cols), "), some groups will be colourless",
                  sep = ""))
  }
  output <- cols[as.integer(group)]
  legend <- table(output)[cols]
  legend[is.na(legend)] <- 0
  if (get) {
    expanded$col <- output
    attr(expanded, "legend") <- paste(legend, "/", sum(legend),
                                      sep = "")
    attr(expanded, "fill") <- cols
    return(expanded)
  }
  else {
    attr(output, "legend") <- paste(legend, "/", sum(legend),
                                    sep = "")
    attr(output, "fill") <- cols
    return(output)
  }
}

#' @rdname colourByValueMltp
#' @export
#'
###############################################################################
#colourByConservationMltp
###############################################################################
colourByConservationMltp <- function(helix,msa,cols,get=FALSE,top.name = "FALSE",sort = TRUE,comp.seq = FALSE){
  if(!is.helix.mltp(helix)){
    stop("invalid helix input")
  }
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }
  #restore warnings
  #options(warn = -1)

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

  output <- colourBy.internal(helixConservation(helix.fake[,1:4], msa.prep.equal), 1, 0, cols)
  if (get) {
    helix$col <- output
    attr(helix, "legend") <- attr(output, "legend")
    attr(helix, "fill") <- attr(output, "fill")
    return(helix)
  } else {
    return(output)
  }
}


#' @rdname colourByValueMltp
#' @export
#'
###############################################################################
#colourByCanonicalMltp
###############################################################################
colourByCanonicalMltp <- function(helix,msa,cols,get=FALSE,top.name = "FALSE",sort = TRUE,comp.seq = FALSE){
  if(!is.helix.mltp(helix)){
    stop("invalid helix input")
  }
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }
  #restore warnings
  #options(warn = -1)

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

  output <- colourBy.internal(helixCanonical(helix.fake[,1:4], msa.prep.equal), 1, 0, cols)
  if (get) {
    helix$col <- output
    attr(helix, "legend") <- attr(output, "legend")
    attr(helix, "fill") <- attr(output, "fill")
    return(helix)
  } else {
    return(output)
  }
}


#' @rdname colourByValueMltp
#' @export
#'
###############################################################################
#colourByCountMltp
###############################################################################
colourByCountMltp <- function(helix,cols,counts, get = FALSE){
  if (!is.helix.mltp(helix)) {
    stop("Not a valid helix data frame, aborting")
  }
  if (missing(cols)) {
    cols <- defaultPalette()
  }
  if (missing(counts)) {
    output <- as.character(cut(seq_len(nrow(helix)),
                               breaks = seq(1,nrow(helix), length.out = length(cols) + 1),
                               include.lowest = TRUE,labels = cols))
  }else {
    if (length(cols) != length(counts)) {
      stop(paste("Length of counts (", length(counts),
                 ") does not match length of cols (", length(cols),
                 ")", sep = ""))
    }
    output <- as.character(rep(cols, counts))
  }

  legend <- table(output)[cols]
  legend[is.na(legend)] <- 0
  if (get) {
    helix$col <- NA
    n <- min(length(output), nrow(helix))
    helix$col[1:n] <- output[1:n]
    attr(helix, "legend") <- paste(legend, "/", sum(legend),
                                   sep = "")
    attr(helix, "fill") <- cols
    return(helix)
  }
  else {
    attr(output, "legend") <- paste(legend, "/", sum(legend),
                                    sep = "")
    attr(output, "fill") <- cols
    return(output)
  }
}


#' @rdname colourByValueMltp
#' @export
#'
###############################################################################
#colourByCovariationMltp
###############################################################################
colourByCovariationMltp <- function(helix,msa,cols,get=FALSE,top.name = "FALSE",sort = TRUE,comp.seq = FALSE){
  if(!is.helix.mltp(helix)){
    stop("invalid helix input")
  }
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }
  #restore warnings
  #options(warn = -1)

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

  output <- colourBy.internal(helixCovariation(helix.fake[,1:4], msa.prep.equal), 2, -2, cols)
  if (get) {
    helix$col <- output
    attr(helix, "legend") <- attr(output, "legend")
    attr(helix, "fill") <- attr(output, "fill")
    return(helix)
  } else {
    return(output)
  }
}



###############################################################################
#colourByBasepairFrequency
###############################################################################
basepairFrequencyMltp <- function(helix) {

  helix <- copy(helix)

  if (!is.helix.mltp(helix)) {
    stop("Not a valid helix data frame, aborting")
  }
  if (nrow(helix) == 0) {
    return(helix)
  }else {
    basepairs <- expandHelixMltp(helix)[, c("i", "j")]
  }
  counts <- table(paste(basepairs$i, basepairs$j))
  pos <- as.integer(unlist(strsplit(names(counts), " ")))
  odds <- seq(1, length(pos), by = 2)
  output <- data.table(i = pos[odds], j = pos[odds + 1], length = 1,
                       value = as.integer(counts),ent.i= expandHelixMltp(helix)[,"ent.i"],
                       ent.j= expandHelixMltp(helix)[,"ent.j"])

  output <- output[order(-output$value, output$i, output$j)]
  output <- as.helix.mltp(output, attr(helix, "length"))

  return(output)
}


#' @rdname colourByValueMltp
#' @export
#'
###############################################################################
#colourByBasepairFrequencyMltp
###############################################################################
colourByBasepairFrequencyMltp <- function(helix, cols, get = TRUE){

  if (!is.helix.mltp(helix)) {
    stop("Not a valid helix data frame, aborting")
  }
  if (missing(cols)) {
    cols <- defaultPalette()
  }

  freq <- colourByValueMltp(basepairFrequencyMltp(helix), cols, get = TRUE)

  if (get) {
    return(freq)
  }
  else {
    output <- freq$col
    attr(output, "legend") <- attr(freq, "legend")
    attr(output, "fill") <- attr(output, "fill")
    return(output)
  }

}





