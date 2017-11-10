################################################################################
# convert meters to lat/lon degrees and vice versa
#
# reference: https://knowledge.safe.com/articles/725/calculating-accurate-length-in-meters-for-lat-long.html

#### degrees to meters #########################################################
d2m <- function(degrees, reflat=45.0, type='avg') {
  mpd <- .mpd(reflat) 
  if(type=='lon') {
    return(degrees * mpd[1])
  }
  if(type=='lat') {
    return(degrees * mpd[2])
  }
  if(type=='avg') {
    return(degrees * mean(mpd))
  }
}

#### meters to degrees #########################################################
m2d <- function(meters, reflat=45.0, type='avg') {
  mpd <- .mpd(reflat)
  if(type=='lon') {
    return(meters / mpd[1])
  }
  if(type=='lat') {
    return(meters / mpd[2])
  }
  if(type=='avg') {
    return(meters / mean(mpd))
  }
}

#### calculate (average) cell area #############################################
getAreaMeters <- function(aRaster) {
  rs <- res(aRaster)
  centerLatitude <- mean(extent(aRaster)[3:4]) 
  if(isLonLat(aRaster)) {
    xlen <- d2m(rs[1], reflat=centerLatitude, type='lon')
    ylen <- d2m(rs[2], reflat=centerLatitude, type='lat')
  } else {
    xlen <- rs[1]
    ylen <- rs[2]
  }
  return(xlen*ylen)
}

# service function to calculate 'meters per degree' ############################
.mpd <- function(reflat) {
  # convert reference latitude into radians
  rlat = reflat*pi/180
  # get how many meters we have in a latitude degree at reference latitude
  latMetersPerDegree <- 111132.92 - 559.82 * cos(2 * rlat) + 1.175 * cos(4 * rlat)
  # get how many meters we have in a longitude degree at reference latitude
  lonMetersPerDegree <- 111412.84 * cos(rlat) - 93.5 * cos(3 * rlat)
  #cat("Converting", meters, "m at", reflat, "degrees latitude:\n")
  #cat("\tlon:", lonMetersPerDegree, " - lat:", latMetersPerDegree, "\n")
  return(c(lonMetersPerDegree, latMetersPerDegree))
}
## EOF ###