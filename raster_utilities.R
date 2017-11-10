###############################################################################
# UAGRA R Scripts                                          raster_utilities.R
###############################################################################
# MAXENT.R - Commodity functions to prepare raster dataset for habitat 
#            selection analyses (e.g. MAXENT, home range analysis)
###############################################################################
#
# Version 1.5
#
# Description: functions to deal with raster data
#
# Usage: in your script, just source() this file: 
#        See function bodies for details.
#
# Requires: rgdal, raster
#
# Returns: 
#
###############################################################################
# created prea@uninsubria.it 20410210
# updated prea@uninsubria.it 20170127
# revision history:
#   prea 20170127 added parseLandsatName function
#   prea 20170317 happy birthday! modified makePercentRaster to handle categorical rasters
#   prea 20170320 modified parseLandsatName to have also month and day
#   prea 20170324 1.5 possibily, check whether r.mess worksadded MESS and MoD functions
#   <author name> <timestamp> - <concise description of the changes made, one per line>
#
###############################################################################
# 
#    Copyright (C) 2014-2017 Damiano G. Preatoni
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    A copy of the GNU General Public License is available at 
#    http://www.gnu.org/licenses/gpl.txt, or can be requested to the Free
#    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#    Boston, MA 02110-1301 USA.
#
###############################################################################

library(raster) # also includes sp
library(rgeos)


################################################################################
# function to make a circular weights matrix of given radius and resolution
# NB radius must me an _even_ multiple of res!
# http://scrogster.wordpress.com/2012/10/05/applying-a-circular-moving-window-filter-to-raster-data-in-r/
make_circ_filter <- function(radius, res) {
  circ_filter <- matrix(NA, nrow=1+(2*radius/res), ncol=1+(2*radius/res))
  dimnames(circ_filter)[[1]] <- seq(-radius, radius, by=res)
  dimnames(circ_filter)[[2]] <- seq(-radius, radius, by=res)
  sweeper <- function(mat) {
    for(row in 1:nrow(mat)){
      for(col in 1:ncol(mat)){
        dist <- sqrt((as.numeric(dimnames(mat)[[1]])[row])^2 + (as.numeric(dimnames(mat)[[1]])[col])^2)
        if(dist<=radius) { mat[row, col] <- 1 }
      }
    }
    return(mat)
  }
  out<-sweeper(circ_filter)
  return(out)
}


################################################################################
#### creates a frequency map for each category starting from a categorical raster
makePercentRasters <- function(aRaster, radius, nodata=0, parallel=FALSE, verbose=FALSE) {
  ptm <- proc.time()
  if(!inherits(aRaster, "Raster")) {
    stop("Input raster must be a 'Raster' object. Aborting.")
  }
  # make a list of raster values, assuming it is a categorical raster
  if(is.factor(aRaster)) {
    cats <- levels(aRaster)[[1]]
    categories <- cats$ID
    labels <- cats$levels
  } else { # brute force approach for non-explicitly categorical rasters
    categories <- unique(aRaster)
    categories <- categories[!categories==nodata]
    labels <- categories
  }
  # get raster spetial resolution
  resolution <- res(aRaster)
  if(!all.equal(resolution[1], resolution[2])) {
    stop("Input raster must have square cells. Aborting.")
  }
  resolution <- resolution[1]
  # calculate moving window size in cells (must be odd!), minimum size is 3x3
  wsize <- as.integer(ceiling(radius / resolution))
  if(wsize%%2==0) wsize <- wsize+1
  if(wsize < 3) wsize <- 3
  radius <- ((resolution * wsize)/2) + (resolution/2)
  filter <- make_circ_filter(radius, resolution)
  filterSum <- sum(filter, na.rm=TRUE)
  if(verbose) {cat('moving window size is', wsize, 'cells, number of cells in circular filter:', filterSum, '\n')}
  outStack <- stack()
  calcPercent <- function(map, cat, filter) {
    idxRast <- calc(map, function(x) x==cat)
    pctRast <- focal(idxRast, filter, fun=function(x) sum(x, na.rm=TRUE)/sum(filter, na.rm=TRUE))
    return(pctRast)
  }
  if(parallel) {
    if(verbose) {cat('Going parallel...\n')}
    require(parallel)
    jobs <- lapply(categories, function(x) mcparallel(calcPercent(aRaster, x, filter), name=x)) 
    outStack <- mccollect(jobs)
    outStack <- stack(outStack)
  } else {
    for(c in categories) {
      if(verbose) {cat('Calculating percentage map for category', c, '\n')}
      idxRaster <- calc(aRaster, fun=function(x) x==c)
      pctRaster <- focal(idxRaster, filter, fun=function(x) sum(x, na.rm=TRUE)/filterSum)
      outStack <- addLayer(outStack, pctRaster)
    }
  }
  names(outStack) <- as.character(labels)
  if(verbose) {print(proc.time() - ptm)}
  invisible(outStack)
} 

################################################################################
#### gets a raster, crops to a polygon and reprojects in the same CRS as the polygon
cropAndReproject <- function(aCover, aCutter, method="bilinear") {
  stopifnot(inherits(aCover, "Raster") | inherits(aCover, "SpatialPolygons"))
  stopifnot(inherits(aCutter, "SpatialPolygons"))    
  # get cutter CRS
  outCRS <- CRS(projection(aCutter))
  covCRS <- CRS(projection(aCover))
  if(projection(aCutter)!=projection(aCover)) { # reproject cutter
    tmpCutter <- spTransform(aCutter, covCRS)
  } else {
    tmpCutter <- aCutter
  }
  # crop and reproject
  if(inherits(aCover, "Raster")) { # we have to process a raster
    cr <- crop(aCover, tmpCutter)
    r <- projectRaster(cr, crs=outCRS, method=method)
  }
  if(inherits(aCover, "SpatialPolygons")) { # process a vector
    cr <- gIntersection(aCover, tmpCutter, byid=TRUE, drop_lower_td=TRUE)
    r <- spTransform(cr, crs=outCRS)
  }
  invisible(r)
}

################################################################################
#### 'reshapes' a raster according to a 'template' raster, matching projection, 
#### extent and cell size
alignRaster <- function(aRaster, aTemplate, method="bilinear", verbose=FALSE) {
  stopifnot(inherits(aRaster, "Raster"))
  stopifnot(inherits(aTemplate, "Raster"))
  # projection
  if(verbose) cat("checking projection...")
  if(projection(aRaster)!=projection(aTemplate)) {
    if(verbose) cat(" reprojecting")
    aRaster <- projectRaster(aRaster, aTemplate, method=method)
  }
  if(verbose) cat(" Done.\n")
  # extent
  if(verbose) cat("matching extent")
  if(extent(aRaster)>extent(aTemplate)) {
    if(verbose) cat(" cropping...")
    aRaster <- crop(aRaster, aTemplate)
  } else {
    if(verbose) cat(" extending...")
    aRaster <- extend(aRaster, aTemplate)
  }
  if(verbose) cat(" Done.\n")
  # cell size
  if(verbose) cat("matching spatial resolution")
  aRaster <- resample(aRaster, aTemplate, method=method)
  if(verbose) cat(" Done.\n")
  invisible(aRaster)
}

parseLandsatName <- function(filename) {
  stopifnot(is.character(filename))
  # as per https://landsat.usgs.gov/what-are-naming-conventions-landsat-scene-identifiers
  # LXSPPPRRRYYYYDDDGSIVV
  # where:
  #  L - Landsat (constant)
  #  X - Sensor (C: OLI/TIRS; O: OLI, T: TIRS; E: ETM; T: TM, M: MSS)
  #  S - Satellite (07: Landsat 7; 08: Landsat 8)
  #  PPP - Path
  #  RRR - Row
  #  YYYY - Year
  #  DDD - Julian day
  #  GSI - Ground Station ID
  #  VV - Archive version number
  sensor <- substr(filename,2,2)
  satellite <- substr(filename,3,3)
  path <- substr(filename,4,6)  
  row <- substr(filename,7,9)  
  year <- substr(filename,10,13)   
  julday <- substr(filename,14,16)   
  theDate <- strptime(paste(year, julday), "%Y %j")
  month <- as.character(format(theDate, format="%m"))
  day <- as.character(format(theDate, format="%d"))
  date <- as.character(theDate)
  ground.station <- substr(filename,17,19)   
  archive <-substr(filename,20,22)   
  result <- c(sensor, satellite, path, row, year, month, day, julday, date, ground.station, archive)
  names(result) <- c('sensor', 'satellite', 'path', 'row', 'year', 'month', 'day', 'julday', 'date', 'ground.station', 'archive')
  invisible(result)
}

# possibly ot nttede, if r.mess in GRASS7 works.
#### calculates MESS in the same way as dismo::mess does
# credits to Jean-Pierre Rossi https://rossijeanpierre.wordpress.com/2012/08/13/computing-the-multivariate-environmental-similarity-surfaces-mess-index-in-r/
# MESS <- function(X, V, full=TRUE) {
#   messi<-function(p, v) {
#     niinf<-length(which(v<=p))
#     f    if(f==0) simi<-100*(p-min(v))/(max(v)-min(v))
#     if(0    if(50<=f & f<100) simi<-2*(100-f)
#        if(f==100) simi<-100*(max(v)-p)/(max(v)-min(v))
#        return(simi)
#   }
#   E<-extract(x=X,y=1:ncell(X))
#   r_mess<-X
#   for (i in 1:(dim(E)[2])) {
#     e<-data.frame(E[,i]) ; v<-V[,i]
#     r_mess[[i]][]<-apply(X=e, MARGIN=1, FUN=messi, v=v)
#   }
#   rmess<-r_mess[[1]]
#   E<-extract(x=r_mess,y=1:ncell(r_mess[[1]]))
#   rmess[]<-apply(X=E, MARGIN=1, FUN=min)
#   if(full==TRUE) {
#     out     layerNames(out)<-c(layerNames(X),"mess")
#   }
#   if(full==FALSE) out return(out)
# }
# 

