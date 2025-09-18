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
#readHelixMltp
###############################################################################
#' Read Helix File for Multiple entities
#'
#' @description Reads Helix file from source file into data.table format
#'
#' @param file file name which has helix in multiple entities
#' @author Volodymyr Tsybulskyi
#' @usage readHelixMltp(file)
#' @details
#' \strong{Helix:} Files start with a header line beginning with # followed by the sequence length,
#' followed by a six-column tab-delimited table (with column names), where each row
#' corresponds to a helix in the structure. The four columns are i and j for the
#' left-most and right-most basepair positions respectively, the length of the
#' helix (converging inwards from i and j, an arbitrary value assigned to the helix,
#' ent.i and ent.j entity for left-most and right-most basepair positions respectively.
#'
#' @examples
#'
#'file <- system.file("extdata", "known.helix.mltp.txt", package = "R4RNA")
#'helix <- readHelixMltp(file)
#'head(helix)
#'
#' @export
#'
readHelixMltp <- function(file){
  #open connection for line with length extraction
  helix_file <- file(file,"r")
  len_group <- strsplit(sub("#(\\s+)", "",substring(readLines(helix_file,n=1),
                                                    2)),"\\s+")[[1]]
  close(helix_file)

  if(length(len_group) == 0){
    stop("Could not read valid length, check helix file")
  }
  #split by delimeter ":"
  tmp <- sapply(len_group, function(x) strsplit(x,split = "[:]")[[1]],
                USE.NAMES = FALSE)
  if(nrow(tmp)!=2){
    stop("Names of entities has : inside. Please remove them")
  }
  #creating dta frame for future attribute
  group <- tmp[1,]
  width <- as.integer(tmp[2,])
  attribute <- data.table(group,width)
  #Start to read main table
  helix <- fread(file)
  changeCols <- c("length")
  helix <- helix[,(changeCols):= lapply(.SD, as.integer), .SDcols = changeCols]
  changeCols <- c("i","j","value")
  helix <- helix[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols]
  #test and flip i and j if it will be needed
  helix_replace <- copy(helix[ent.i==ent.j & i>j])
  if (helix_replace[,.N]!=0) {
    warning("\nFlip coordinates i and j for i>j in same entity")
    #flip positions
    tmp_i <- helix_replace[,i]
    tmp_j <- helix_replace[,j]
    helix_replace[,i := tmp_j]
    helix_replace[,j := tmp_i]
    #find indeces and replace them with flipped data.table
    indeces_replace <- helix[,.I[ent.i==ent.j & i>j]]
    helix[c(indeces_replace)] <- helix_replace
  }
  #check if length is higher than it supposed to be by length line
  len_i <- as.vector(sapply(helix[,ent.i], function(x)
    attribute[group==x,as.numeric(width)]))
  len_j <- as.vector(sapply(helix[,ent.j], function(x)
    attribute[group==x,as.numeric(width)]))
  helix[,len.i:= len_i]
  helix[,len.j:= len_j]
  if(length(helix[,.N[i>len.i | j>len.j]]) != 0){
    stop("Interactions out of bounds\nPlease check helix file")
  }
  #attach attribute to helix file
  attr(helix,"length") <- attribute
  return(helix[,1:6])
}



###############################################################################
#is.helix.mltp
###############################################################################
#' Check a Helix Data Table
#'
#'@description
#'Functions to check whether a structure is a valid helix data table and generate helix from data.table.
#'@author Volodymyr Tsybulskyi
#'@param helix helix structure to check
#'
#' @export
#'
is.helix.mltp <- function(helix){
  if("col" %in% colnames(helix)){
    helix.to.test<- copy(helix[,!"col"])
  }else{
    helix.to.test <- copy(helix)
  }

  valid <- TRUE

  if(!is.data.table(helix.to.test)){
    warning("Not a valid data table")
    valid <- FALSE
  }
  if(ncol(helix.to.test)!=6){
    warning("Should have 6 columns")
    valid <- FALSE
  }
  if(!"length" %in% colnames(helix.to.test)) {
    warning("No length attribute")
    valid <- FALSE
  }
  if (!all(lapply(helix.to.test, class)[1:6] == c("numeric", "numeric",
                                                  "integer", "numeric",
                                                  "character","character")) &
      !all(lapply(helix.to.test, class)[1:6] == c("integer", "integer",
                                                  "integer", "numeric",
                                                  "character","character"))) {
    warning("Columns of invalid class")
    valid <- FALSE
  }
  if (!all(colnames(helix.to.test)[1:6] == c("i", "j", "length",
                                             "value","ent.i","ent.j"))) {
    warning("Invalid column names")
    valid <- FALSE
  }
  if (length(attr(helix.to.test,"length"))==0){
    warning("no attribute attachment")
    valid <- FALSE
  }
  return(valid)
}


###############################################################################
#as.helix.mltp
###############################################################################
#' @rdname is.helix.mltp
#' @export
#'
as.helix.mltp <- function(x, length = NULL){
  x <- as.data.table(x)
  if (nrow(x) == 0) {
    x <- data.table(i = vector(mode = "numeric"), j = vector(mode = "numeric"),
                    length = vector(mode = "integer"), value = vector(mode = "numeric"),
                    ent.i = vector(mode = "character"),ent.j = vector(mode = "character"))
  }
  if (ncol(x) < 6) {
    stop("Expecting data.frame-like structure with at least 6 columns")
  }
  names <- colnames(x)
  names[1:6] <- c("i", "j", "length", "value","ent.i","ent.j")
  colnames(x) <- names
  rownames(x) <- NULL

  x$i <- as.numeric(x$i)
  x$j <- as.numeric(x$j)
  x$length <- as.integer(x$length)
  x$value <- as.numeric(x$value)
  x$ent.i <- as.character(x$ent.i)
  x$ent.j <- as.character(x$ent.j)

  if (is.null(length)) {
    groups <- unique(c(x$ent.i,x$ent.j))
    lengths <- c()
    for (k in 1:length(groups)) {
      lengths <- c(lengths,as.integer(max(c(x[ent.i==groups[k],i],x[ent.j==groups[k],j]))))
    }
    length <- data.table(group = groups,width = lengths)
  }

  attr(x, "length") <- length
  return(x)
}


###############################################################################
#function to expand row of helix
###############################################################################
small.expand <- function(helix.row){
  if(nrow(helix.row)>1){
    stop("invalid input for function small.expand")
  }
  #multiply row to length
  helix.small.exp <- helix.row[rep(1,length[1])]
  #make calculation for interactions
  helix.small.exp[,i := i+(.I-1)]
  helix.small.exp[,j := j-(.I-1)]
  #replace length by 1
  helix.small.exp[,length:= rep(as.integer(1),helix.small.exp[,.N])]
  return(helix.small.exp)
}


###############################################################################
#expand helix
###############################################################################
#' @export
#'
expandHelixMltp <- function(helix){
  if (!is.helix.mltp(helix)) {
    stop("Input is not a valid helix, aborting")
  }else{
    if (nrow(helix) == 0 || all(helix$length == 1)) {
      return(helix)
    }
  }
  if(helix[ent.i!=ent.j & length>1,.N] > 0){
    stop("You provide trans interactions with length more than 1")
  }
  helix.input <- copy(helix)
  len <- attr(helix,"length")
  #define rows to expand
  helix.to.exp <- helix.input[length>1]
  #delete all rows from input data.table
  helix.out <- copy(helix.input[length==1])
  #expand remain helices
  helix.exp <- helix.to.exp[0,]
  #expand rows
  for (i in 1:helix.to.exp[,.N]) {
    helix.exp <- rbind(helix.exp,small.expand(helix.to.exp[i]))
  }
  #bind to main body
  helix.out <- rbind(helix.out,helix.exp)
  attr(helix.out,"length") <- attr(helix,"length")
  return(helix.out)
}



###############################################################################
#collapse helix
###############################################################################
#' @export
#'
collapseHelixMltp <- function(helix, number = FALSE){
  if (!is.helix.mltp(helix)) {
    stop("Input is not a valid helix, aborting")
  }else {
    if (nrow(helix) == 0) {
      return(helix)
    }
    helix <- expandHelixMltp(helix)
    helix <- helix[!duplicated(helix), ]
  }

  length <- attr(helix, "length")
  tmp <- NA
  if (any(is.na(helix$value))) {
    tmp <- max(helix$value, na.rm = TRUE) + 1
    helix$value[is.na(helix$value)] <- tmp
  }

  #helix to work with
  tmp_exp <- helix[ent.i==ent.j]

  sums <- rowSums(tmp_exp[, 1:2])
  if (number) {
    order <- order(tmp_exp$number, sums, tmp_exp$i)
  }else {
    order <- order(tmp_exp$value, sums, tmp_exp$i)
  }

  tmp_exp <- tmp_exp[order, ]
  sums <- sums[order]
  vend <- rep(FALSE, nrow(tmp_exp))
  vsta <- rep(FALSE, nrow(tmp_exp))
  if (number) {
    ends <- findInterval(unique(tmp_exp$number), tmp_exp$number)
    vend[ends] <- TRUE
    vsta[c(0, ends[-length(ends)]) + 1] <- TRUE
  }else {
    ends <- findInterval(unique(tmp_exp$value), tmp_exp$value)
    vend[ends] <- TRUE
    vsta[c(0, ends[-length(ends)]) + 1] <- TRUE
  }

  ends <- cumsum(rle(sums)$lengths)
  vend[ends] <- TRUE
  vsta[c(0, ends[-length(ends)]) + 1] <- TRUE
  ends <- which(c(diff(tmp_exp$i), 0) != 1)
  vend[ends] <- TRUE
  vsta[c(0, ends[-length(ends)]) + 1] <- TRUE
  starts <- which(vsta)
  ends <- which(vend)
  counts <- which(vend) - which(vsta) + 1
  ungap <- (tmp_exp$i[vsta] == tmp_exp$i[vend] - counts + 1) &
    (tmp_exp$j[vsta] == tmp_exp$j[vend] + counts - 1)
  if (length(which(ungap == FALSE)) > 0) {
    stop("This error should never occur... inform developer if it does")
  }

  output <- tmp_exp[vsta, ]
  output$length <- counts

  output <- rbind(output, copy(helix[ent.i!=ent.j | ent.j!=ent.i]))

  return(as.helix.mltp(output, length))

}


###############################################################################
#writeHelixMltp
###############################################################################
#' @export
#'
writeHelixMltp <- function(file,helix,uniq){
  length <- paste(uniq$group,uniq$width,sep = ":")
  length[1] <- paste("#",length[1],sep = "")
  length <- paste(length,collapse = " ")
  names <- paste(c("i", "j","length", "value", "ent.i", "ent.j"), collapse = "\t")

  write(length,file=file)
  write(names,file=file,append=TRUE)
  fwrite(file = file,x = helix,append = TRUE,sep = "\t")
}


###############################################################################
#StrawToHelix
###############################################################################
#' @export
#'
StrawToHelix <- function(file,chr1,chr2,scale){
  #read file
  obj1 <- read.table(file, header=FALSE)
  colnames(obj1) <- c("i","j","value")

  obj1 <- as.data.table(obj1)
  #scale object for better representation
  obj1.part1 <- obj1[i >0][j>0]
  obj1.part1 <- obj1.part1[,i:= i / scale][,j:= j / scale]

  obj1.part2 <- obj1[i==0][j>0]
  obj1.part2 <- obj1.part2[,j:= j / scale]

  obj1.part3 <- obj1[i>0][j==0]
  obj1.part3 <- obj1.part3[,i:= i / scale]

  obj1.part4 <- obj1[i==0][j==0]

  obj1 <- rbind(obj1.part4,obj1.part3,obj1.part2,obj1.part1)

  obj1[,length := as.integer(1)][,ent.i := chr1][,ent.j := chr2]
  setcolorder(obj1, c("i", "j", "length", "value", "ent.i", "ent.j"))
  obj1 <- obj1[,value := as.numeric(value)]
  obj1 <- obj1[,i := i+1][,j:=j+1]
  obj1 <- obj1[order(value)]

  uniq <- data.table(group = c(chr1,chr2),width = c(max(obj1$i),max(obj1$j)))
  attr(obj1,"length") <- uniq

  return(obj1)
}




