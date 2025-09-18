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


#'Compute statistics for a multiple sequence alignments
#'
#'@description
#'Functions to compute covariation, percent identity conservation,
#'and percent canonical basepairs given a multiple sequence alignment
#'and optionally a secondary structure. Statistics can be computed for a
#'single base, basepair, helix or entire alignment.
#'
#' @param msa Multiple sequence alignment as an array of named characters,
#' all of equal length. Typically output of readBstringSet.
#' @param helix Helix data.tables, with the 6 mandatory columns. "col" refer to a styling column, and will be used for styling the helix. See example for styling usage.
#' @param top.name Name of entity, which will be set as top name and will refer to first entity in plot
#' @param sort sort entities by alphabetic order
#' @param comp.seq Option to compare sequences from different entities by their "name" and provide
#' proper covariance colouring in future
#'
#' @details
#' Conservation values have a range of [0, 1], where 0 is the absence of
#' primary sequence conservation (all bases different), and 1 is full primary
#' sequence conservation (all bases identical).
#'
#' Canonical values have a range of [0, 1], where 0 is a complete lack of
#' basepair potential, and 1 indicates that all basepairs are valid
#'
#' Covariation values have a range of [-2, 2], where -2 is a complete
#' lack of basepair potential and sequence conservation, 0 is complete
#' sequence conservation regardless of basepairing potential, and 2 is a
#' complete lack of sequence conservation but maintaining full basepair potential.
#'
#' helix values are average of base/basepair values, and the alignment
#' values are averages of helices or all columns depending on whether
#' the helix argument is required.
#'
#' @export
#'
###############################################################################
#alignmentCanonicalMltp
###############################################################################
alignmentCanonicalMltp <- function(msa,helix,top.name = "FALSE",sort = TRUE,comp.seq = FALSE){

  if(!is.helix.mltp(helix)){
    stop("invalid helix input")
  }
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }

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

  helix.fake <- helix.fake[,1:4]

  percentCanonical <- 0
  for (bp.num in 1:nrow(helix.fake)) {
    bp.i <- helix.fake[bp.num, "i"]
    bp.j <- helix.fake[bp.num, "j"]
    percentCanonical <- percentCanonical + basepairCanonical(msa,
                                                             bp.i, bp.j)
  }

  return(percentCanonical/nrow(helix.fake))
}


#' @rdname alignmentCanonicalMltp
#' @export
#'
###############################################################################
#alignmentCovariationMltp
###############################################################################
alignmentCovariationMltp <- function(msa,helix,top.name = "FALSE",sort = TRUE,comp.seq = FALSE){

  if(!is.helix.mltp(helix)){
    stop("invalid helix input")
  }
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }

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

  helix.fake <- helix.fake[,1:4]

  covariation <- sum(apply(helix.fake, 1, function(bp, msa.prep.equal) {
    basepairCovariation(msa, bp["i"], bp["j"])
  }, msa))
  return(covariation/nrow(helix))

}

#' @rdname alignmentCanonicalMltp
#' @export
#'
###############################################################################
#helixCovariationMltp
###############################################################################
helixCovariationMltp <- function(msa,helix,top.name = "FALSE",sort = TRUE,comp.seq = FALSE){

  if(!is.helix.mltp(helix)){
    stop("invalid helix input")
  }
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }

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

  helix.fake <- helix.fake[,1:4]

  ex <- copy(expandHelix(helix.fake))
  num <- rep(1:nrow(helix.fake), times = helix.fake$length)
  cov <- c()
  for (i in 1:nrow(ex)) {
    cov[i] <- basepairCovariation(msa.prep.equal, ex$i[i], ex$j[i])
  }
  pt <- 1
  out = c()
  for (i in 1:nrow(helix.fake)) {
    out[i] <- mean(cov[seq(pt, length.out = helix.fake$length[i])])
    pt <- pt + helix.fake$length[i]
  }
  return(out)

}




#' @rdname alignmentCanonicalMltp
#' @export
#'
###############################################################################
#helixConservationMltp
###############################################################################
helixConservationMltp <- function(msa,helix,top.name = "FALSE",sort = TRUE,comp.seq = FALSE){

  if(!is.helix.mltp(helix)){
    stop("invalid helix input")
  }
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }

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

  helix.fake <- helix.fake[,1:4]

  ex <- copy(expandHelix(helix.fake))
  num <- rep(1:nrow(helix.fake), times = helix.fake$length)
  cons <- c()
  for (i in 1:nrow(ex)) {
    cons[i] <- basepairConservation(msa.prep.equal, ex$i[i], ex$j[i])
  }
  pt <- 1
  out = c()
  for (i in 1:nrow(helix.fake)) {
    out[i] <- mean(cons[seq(pt, length.out = helix.fake$length[i])])
    pt <- pt + helix.fake$length[i]
  }
  return(out)


}


#' @rdname alignmentCanonicalMltp
#' @export
#'
###############################################################################
#helixCanonicalMltp
###############################################################################
helixCanonicalMltp <- function(msa,helix,top.name = "FALSE",sort = TRUE,comp.seq = FALSE){

  if(!is.helix.mltp(helix)){
    stop("invalid helix input")
  }
  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }

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

  helix.fake <- helix.fake[,1:4]

  ex <- copy(expandHelix(helix.fake))
  num <- rep(1:nrow(helix.fake), times = helix.fake$length)
  cano <- c()
  for (i in 1:nrow(ex)) {
    cano[i] <- basepairCanonical(msa.prep.equal, ex$i[i], ex$j[i])
  }
  pt <- 1
  out = c()
  for (i in 1:nrow(helix.fake)) {
    out[i] <- mean(cano[seq(pt, length.out = helix.fake$length[i])])
    pt <- pt + helix.fake$length[i]
  }
  return(out)

}


#' @rdname alignmentCanonicalMltp
#' @export
#'
###############################################################################
#baseConservationMltp
###############################################################################
baseConservationMltp <- function(msa,pos,GroupName){

  if (is(msa, "XStringSet")) {
    msa <- as.character(msa)
  }

  #####Extract sequences and make msa as 1 BStringSet
  #####For next using on default functions
  #transfor fasta to data.table format
  msa.data.table <- bstring.to.data.table(msa,comp.seq = comp.seq)
  groups <- unique(msa.data.table$group)

  if (missing(GroupName)) {
    GroupName <- groups
  }
  if (!GroupName %in% groups) {
    stop("Your GroupName not in msa names\n")
  }

  out = data.table(group = vector(mode = "character"),conservation = vector(mode = "numeric"))

  for (k in 1:length(GroupName)) {
    nameGroup <- GroupName[k]
    x = msa.data.table[group %like% nameGroup]
    if (pos > unique(x$width)) {
      stop("Position out of bounds.")
    }
    msa <- BStringSet(x = x$sequence)
    value <- columnPercentIdentity(substr(msa, pos, pos))

    out <- rbind(out,data.table(group = nameGroup,conservation = value))
  }

  return(out)
}



#'Logical filters of helix by type
#'
#' @description
#'
#' Given a helix data frame, checks if helices are conflicting,
#' duplicating, or overlapping, and returns an array of logicals.
#' See details for exact definition of the three types of events.
#'
#' @details
#' Helices of length greater than 1 are internall expanded into
#' basepairs of length 1, after which the following conditions are evaluated:
#' A conflicting basepair is one where at least one of its two positions is used by either end of another basepair.
#' A duplicating basepair is one where both of its positions are used by both ends of another basepair.
#' An overlapping basepair is one in helix where both of its
#' positions are used by both ends of another basepair in the query structure.
#'
#' In the case of conflicting and duplicating basepairs, for a set of basepairs
#' that satisfies this condition, the basepair situation highest on the data frame
#' will be exempt from the condition. i.e. Say 5 basepairs are all duplicates of each
#' other, the top 1 will return FALSE, while the bottom 4 will return TRUE. This
#' assumes some significant meaning to the ordering of rows prior to using this function.
#' This is to be used with which to filter out basepairs that satisfy these conditions,
#' leaving a set of basepairs free of these events.
#'
#' @param helix A helix data table
#' @param query For \code{isOverlappingHelixMltp}, a helix data structure against
#' which \code{helix} will be checked for overlap against.
#'
#' @export
#'
###############################################################################
#isConflictingHelixMltp
###############################################################################
isConflictingHelixMltp <- function(helix){
  if (!is.helix.mltp(helix)) {
    stop("Not a valid helix data.table, aborting...")
  }
  if (nrow(helix) == 0) {
    return(vector())
  }
  if (nrow(helix) == 1) {
    return(c(FALSE))
  }

  helix <- expandNumberedHelixMltp(helix)
  end <- nrow(helix)
  helix$bad <- 1
  for (i in 1:(end - 1)) {
    if (helix$bad[i] == 1) {
      kill <- which(helix$i[(i + 1):end] == helix$i[i] | helix$j[(i + 1):end] == helix$j[i] |
                      helix$j[(i + 1):end] == helix$i[i] | helix$i[(i + 1):end] == helix$j[i])
      helix$bad[kill + i] <- 0
    }
  }
  return(aggregate(1 - helix$bad, by = list(helix$number), mean)$x)
}


#' @rdname isConflictingHelixMltp
#' @export
#'
###############################################################################
#isDuplicatingHelixMltp
###############################################################################
isDuplicatingHelixMltp <- function(helix){

  if (!is.helix.mltp(helix)) {
    stop("Not a valid helix data.table, aborting...")
  }
  if (nrow(helix) == 0) {
    return(vector())
  }
  if (nrow(helix) == 1) {
    return(c(FALSE))
  }

  expand <- expandNumberedHelixMltp(helix)
  flip <- which(expand$i > expand$j)
  tmp <- expand$i[flip]
  expand$i[flip] <- expand$j[flip]
  expand$j[flip] <- tmp
  expand$duplicate <- duplicated(expand[, c("i", "j","ent.i","ent.j")])
  return(aggregate(expand$duplicate, by = list(expand$number),mean)$x)

}


#' @rdname isConflictingHelixMltp
#' @export
#'
###############################################################################
#isOverlappingHelixMltp
###############################################################################
isOverlappingHelixMltp <- function(helix, query){

  if (!is.helix.mltp(helix) | !is.helix.mltp(query)) {
    stop("One of the inputs was not a valid helix data frame, aborting")
  }
  if (nrow(helix) == 0) {
    return(vector())
  }
  if (nrow(query) == 0) {
    return(rep(0, helix))
  }

  #define which are known and unknown
  setDT(helix)
  setDT(query)
  helix.processed <- helix[, hit := 0][query, hit := 1,
                                       on = .(i,j,ent.i,ent.j)]

  data.return <- helix.processed$hit
  return(data.return)
}
