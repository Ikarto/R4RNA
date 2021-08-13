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






#' Plots helices in diagonal arc diagram
#'
#' @description
#' Plots a helix data frame as an diagonal arc diagram,
#' with styling possible with properly named additional columns on the data frame.
#'
#' @param helix Helix data.tables, with the 6 mandatory columns. "col" refer to a styling column, and will be used for styling the helix. See example for styling usage.
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
#' @param add,append If TRUE, graphical elements are added to the active plot device, else a new plot device is created for the plot.
#' @param arc.lty define lty of arcs
#'
#' @export
#'
###############################################################################
#plotContact
###############################################################################
plotContact <- function(helix,direction = "se",top.name = "FALSE",sort = TRUE,
                        dist.part = .2,stable = FALSE,line = TRUE,
                        arrow = TRUE,col.line = "black",
                        col.arrow = "black",col = "black",
                        shape = "triangle",scale = FALSE,
                        debug = FALSE,add = FALSE,flip = FALSE,arc.lty = 1,append = FALSE,...){

  #reduce data.table warning
  #options(warn = -1)

  if(helix[length>1,.N]!=0){
    stop("Helix is not expanded")
  }
  if(is.null(attr(helix,"length"))){
    stop("No attribute for helix")
  }
  if(!direction %in% c("se","ne","SE","NE")){
    stop("Direction option wrong\n Choose between: 'se','ne'")
  }

  #read len data
  uniq <- attr(helix,"length")
  #sort and define top name
  uniq <- update.uniq.by.group.single(uniq,top.name,sort)

  ##############################################################################
  if(direction == "se" | direction == "SE"){
    #calculate values for uniq
    uniq <- prep.uniq.single(uniq,dist.part)
    full.len <- attr(uniq,"full.length")
    dist.between <- attr(uniq,"dist.between")
    #update helix to new coordinates
    helix <- update.helix.coord.single(helix,uniq)

    diag_full <- seq(max(uniq$end),1)
    for(k in 1:nrow(helix)) {
      helix[k,ent.i.y := diag_full[helix[k,i.coord]]]
      helix[k,ent.j.y := diag_full[helix[k,j.coord]]]
    }
    max.height <- max(uniq$end)
  }
  if(direction == "ne" | direction == "NE"){
    #calculate values for uniq
    uniq <- prep.uniq.single(uniq,dist.part)
    full.len <- attr(uniq,"full.length")
    dist.between <- attr(uniq,"dist.between")
    #update helix to new coordinates
    helix <- update.helix.coord.single(helix,uniq)

    diag_full <- seq(1,max(uniq$end))
    for(k in 1:nrow(helix)) {
      helix[k,ent.i.y := diag_full[helix[k,i.coord]]]
      helix[k,ent.j.y := diag_full[helix[k,j.coord]]]
    }
    max.height <- max(uniq$end)
  }

  if(direction == "nw" | direction == "NW"){
    #calculate values for uniq
    uniq <- prep.uniq.single.reverse(uniq,dist.part)
    full.len <- attr(uniq,"full.length")
    dist.between <- attr(uniq,"dist.between")
    #update helix to new coordinates
    helix <- update.helix.coord.single.reverse(helix,uniq)

    diag_full <- seq(max(c(uniq$end,uniq$start)),1)

    for(k in 1:nrow(helix)) {
      helix[k,ent.i.y := diag_full[helix[k,i.coord]]]
      helix[k,ent.j.y := diag_full[helix[k,j.coord]]]
    }
    max.height <- max(c(uniq$end,uniq$start))
  }

  if(direction == "sw" | direction == "SW"){

    #calculate values for uniq
    uniq <- prep.uniq.single.reverse(uniq,dist.part)
    full.len <- attr(uniq,"full.length")
    dist.between <- attr(uniq,"dist.between")
    #update helix to new coordinates
    helix <- update.helix.coord.single.reverse(helix,uniq)

    diag_full <- seq(1,max(c(uniq$end,uniq$start)))

    for(k in 1:nrow(helix)) {
      helix[k,ent.i.y := diag_full[helix[k,i.coord]]]
      helix[k,ent.j.y := diag_full[helix[k,j.coord]]]
    }
    max.height <- max(c(uniq$end,uniq$start))
  }

  ##############################################################################
  #create space
  if(!add){
    if(append){
      blankPlotMod(full.len+(dist.between*nrow((uniq)))/8,
                   max.height, -1,scale = scale,debug = debug,
                   no.par = TRUE,...)
    }else{
      blankPlotMod(full.len+(dist.between*nrow((uniq)))/8,
                   max.height, -1,scale = scale,debug = debug,...)
    }
  }


  if(direction == "se" | direction == "SE"){
    #plot bases
    for(k in 1:nrow(uniq)) {
      uniq[k,y.start := diag_full[uniq[k,start]]]
      uniq[k,y.end := diag_full[uniq[k,end]]]
    }
    uniq[is.na(uniq)] <- diag_full[1]

    if(line){
      for(i in 1:nrow(uniq)){
        lines(c(uniq[i,start],uniq[i,end]),c(uniq[i,y.start],uniq[i,y.end]),col = col.line)
      }
    }
    if(arrow){
      for (i in 1:uniq[,.N]) {
        x.coord <- c(uniq[i,end]-1,uniq[i,end]+1,uniq[i,end]+dist.between/8)
        y.coord <- c(uniq[i,y.end]-1,uniq[i,y.end]+1,uniq[i,y.end]-dist.between/8)
        polygon(x.coord, y.coord,col = col.arrow)
      }
    }

    plot.arcs.contact(helix = helix,col = col,flip = TRUE,shape = shape)
    plot.arcs.contact(helix = helix,col = col,flip = FALSE,shape = shape)


  }

  if(direction == "ne" | direction == "NE"){
    uniq[start==0,start := 1]
    for(k in 1:nrow(uniq)) {
      uniq[k,y.start := diag_full[uniq[k,start]]]
      uniq[k,y.end := diag_full[uniq[k,end]]]
    }
    uniq[is.na(uniq)] <- diag_full[1]

    if(line){
      for(i in 1:nrow(uniq)){
        lines(c(uniq[i,start],uniq[i,end]),c(uniq[i,y.start],uniq[i,y.end]),col = col.line)
      }
    }
    if(arrow){
      for (i in 1:uniq[,.N]) {
        x.coord <- c(uniq[i,end]-1,uniq[i,end]+1,uniq[i,end]+dist.between/8)
        y.coord <- c(uniq[i,y.end]+1,uniq[i,y.end]-1,uniq[i,y.end]+dist.between/8)
        polygon(x.coord, y.coord,col = col.arrow)
      }
    }


    plot.arcs.contact(helix = helix,col = col,flip = TRUE,shape = shape)
    plot.arcs.contact(helix = helix,col = col,flip = FALSE,shape = shape)

  }

  if(direction == "nw" | direction == "NW"){
    uniq[end==0,end := 1]

    for(k in 1:nrow(uniq)) {
      uniq[k,y.start := diag_full[uniq[k,start]]]
      uniq[k,y.end := diag_full[uniq[k,end]]]
    }
    uniq[is.na(uniq)] <- diag_full[1]
    if(line){
      for(i in 1:nrow(uniq)){
        lines(c(uniq[i,start],uniq[i,end]),c(uniq[i,y.start],uniq[i,y.end]),col = col.line)
      }
    }
    if(arrow){
      for (i in 1:uniq[,.N]) {
        x.coord <- c(uniq[i,end]-1,uniq[i,end]+1,uniq[i,end]-dist.between/8)
        y.coord <- c(uniq[i,y.end]-1,uniq[i,y.end]+1,uniq[i,y.end]+dist.between/8)
        polygon(x.coord, y.coord,col = col.arrow)
      }
    }

    plot.arcs.contact(helix = helix,col = col,flip = TRUE,shape = shape)
    plot.arcs.contact(helix = helix,col = col,flip = FALSE,shape = shape)

  }

  if(direction == "sw" | direction == "SW"){
    uniq[end==0,end := 1]

    for(k in 1:nrow(uniq)) {
      uniq[k,y.start := diag_full[uniq[k,start]]]
      uniq[k,y.end := diag_full[uniq[k,end]]]
    }
    uniq[is.na(uniq)] <- diag_full[1]
    if(line){
      for(i in 1:nrow(uniq)){
        lines(c(uniq[i,start],uniq[i,end]),c(uniq[i,y.start],uniq[i,y.end]),col = col.line)
      }
    }
    if(arrow){
      for (i in 1:uniq[,.N]) {
        x.coord <- c(uniq[i,end]-1,uniq[i,end]+1,uniq[i,end]-dist.between/8)
        y.coord <- c(uniq[i,y.end]+1,uniq[i,y.end]-1,uniq[i,y.end]-dist.between/8)
        polygon(x.coord, y.coord,col = col.arrow)
      }
    }

    plot.arcs.contact(helix = helix,col = col,flip = TRUE,shape = shape)
    plot.arcs.contact(helix = helix,col = col,flip = FALSE,shape = shape)

  }


  #restore warnings
  options(warn = 1)

}
