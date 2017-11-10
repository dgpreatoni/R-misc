################################################################################
# a sample R script to use raster_utilities.R functions.
################################################################################

rm(list=ls())

rootDir <- '/home/gis' # on eira
rootDir <- 'M:'        # on Windows boxes

R_libraries_path <- paste(rootDir,'R_scripts', sep='/')
source(paste(R_libraries_path, 'raster_utilities.R', sep='/'))


# define data directories
# we assume that all data are under a 'rootDir': 
dataDir <- paste(rootDir, 'Scoiattolo/Callosciurus', sep='/') # just to have a raster
BASEDir <- paste(dataDir, 'BASE', sep='/')

vegMapFile <- 'carta_vegetazione_SIT-Fauna/vegx_C4_32632.tif'

# load vegetation map
vegMap <- raster(paste(BASEDir, vegMapFile, sep='/'))

# calculate percentage maps
pct030 <- makePercentRasters(vegMap, radius=30, parallel=FALSE, verbose=TRUE)
pct050 <- makePercentRasters(vegMap, radius=50, parallel=TRUE)
pct100 <- makePercentRasters(vegMap, radius=100, parallel=TRUE)

# raster objects can be exported (e.g. as GeoTIFF) "as usual"
  writeRaster(pct030*100, filename=tempdir(), format="GTiff", bylayer=TRUE, overwrite=TRUE)
