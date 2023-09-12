################################################################################
# preprocessing_EGV.R
#
# Pre-processes raw, partial GIS data to prepare EcoGeographical Variables for
# species distribution modeling.
# EGVs will be stored as a raster::stack in ./BASE/EGV/EGV.Rdata
#
################################################################################
# Version 1.4
# Created: 20181112 prea
# Modified: 20221124 prea
# Version history:
# 1.4 - adapted for stambeccorobie MVA, started using sf instead of sp
# 1.3 - adapted for rischio orso valtellina scenario
# 1.2 - final version, extended AOI
# 1.1 - moved on RStudio server
#
################################################################################

library(raster)
library(rgdal)
library(sf)
library(readxl)
library(parallel)

source(file.path('scripts', 'raster_utilities.R'))

# set up a scratch directrory for rasters
rasterOptions(tmpdir = '/srv/scratch/')

# we want everything in ETRS89/LAEA
theCRS <- CRS("+init=epsg:3035")

# spatial resolution
theResolution <- 100 # meters

# neighborhood radius for 'percentage maps'
theRadius <- 1000 # meters

# name for the EGV stack (will go in ./BASE/EGV)
EGVStackName <- 'EGV.RData'

# reference frame for Area Of Interest
AOI <- st_read(dsn='BASE', layer='AOI 3035') 
AOI <- st_transform(AOI, theCRS)

#### generate an empty 'reference raster'
refRaster <- raster(extent(AOI), resolution=theResolution, crs=theCRS, vals=NA)

#### start up parallel processing to speed up raster calculations
nCores <- parallelly::availableCores(omit=2)
beginCluster(n=nCores) # set up things to run with at least 8 vCPUs, ask BOFH to do that!

########################################################## DEM derived data ####
#### elevation
raw_dem <- raster(file.path('BASE', 'EGV', 'EUDEM_25_Alps.tif'))
# ensure projection is OK(
if(projection(raw_dem) != projection(theCRS)) { # reproject only if needed
  dem <- projectRaster(from=raw_dem, to=refRaster, res=theResolution, method='bilinear')
} else {
  dem <- raw_dem
}
names(dem) <- "elevation"
rm(raw_dem)

#### derive slope, aspect, roughness index
dem_derived <- terrain(dem, opt=c('slope', 'aspect', 'roughness'), unit='degrees', neighbors=8)

#### slope
## calculate slope classes maps (classes as per https://geographyfieldwork.com/SlopeSteepnessIndex.htm)
# 0 to 5 degrees -> level
# 5 to 8.5 degrees -> moderate
# 8.5 to 24 -> strong
# > 24 -> steep
slopeClasses <- matrix(c(0,5,0, 5,8.5,1, 8.5,24,2, 24,90,3), byrow=TRUE, ncol=3)
slope <- ratify(reclassify(dem_derived$slope, rcl=slopeClasses))
slope.lvl <- levels(slope)[[1]]
slope.lvl$slope <- c('level', 'moderate', 'strong', 'steep')
slope.lvl$code <- c(0,1,2,3)
levels(slope) <- slope.lvl
slope <- deratify(slope, 'slope')
slope <- layerize(slope)
names(slope) <- slope.lvl$slope 
#
filter <- make_circ_filter(theRadius, theResolution)
filterSum <- sum(filter, na.rm=TRUE)
slopeStack <- stack()
for(n in names(slope)) {
 pct <- focal(slope[[n]], filter, fun=function(x) sum(x, na.rm=TRUE)/filterSum)
 slopeStack <- addLayer(slopeStack, pct)
}
names(slopeStack) <- paste0('slope-', names(slope))

#### aspect
## calculate aspect classes maps (classes as per https://geographyfieldwork.com/SlopeSteepnessIndex.htm)
# 0 to 22.5 degrees -> N
# 22.5 to 67.5 -> NE
# 67.5 to 112.5 -> E
# 112.5 to 157.5 -> SE
# 157.5 to 202.5 -> S
# 202.5 to 247.5 -> SW
# 247.5 to 292.5 -> W
# 292.5 to 337.5 -> NW
# 337.5 to 360 ->  N (again)
aspectClasses <- matrix(c(0,22.5,1, 22.5,67.5,2, 67.5,112.5,3, 112.5,157.5,4, 157.5,202.5,5, 202.5,247.5,6, 247.5,292.5,7, 292.5,337.5,8, 337.5,360,1), byrow=TRUE, ncol=3)
aspect <- ratify(reclassify(dem_derived$aspect, rcl=aspectClasses))
aspect.lvl <- levels(aspect)[[1]]
aspect.lvl$aspect <- c('flat', 'N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW')
aspect.lvl$code <- c(0,1,2,3,4,5,6,7,8)
levels(aspect) <- aspect.lvl
aspect <- deratify(aspect, 'aspect')
aspect <- layerize(aspect)
names(aspect) <- aspect.lvl$aspect 
#
filter <- make_circ_filter(theRadius, theResolution)
filterSum <- sum(filter, na.rm=TRUE)
aspectStack <- stack()
for(n in names(aspect)) {
  pct <- focal(aspect[[n]], filter, fun=function(x) sum(x, na.rm=TRUE)/filterSum)
  aspectStack <- addLayer(aspectStack, pct)
}
names(aspectStack) <- paste0('aspect-', names(aspect))

#### roughness is already OK

#### pack it all
EGVStack <- stack(dem, aspectStack, slopeStack)

#### cleanup
rm(aspect, aspectStack, slope, slopeStack, dem_derived, dem, dem_derived, aspect.lvl, aspectClasses, slope.lvl, slopeClasses)

################################################### land cover derived data ####
# using a cropped version of 100 m CLC2018_V2018_20b2 raster, obtained by
# gdalwarp -cutline ../AOI\ 3035.shp -crop_to_cutline /srv/archivio/Cartografia/Europa/CORINE_CLC2018/clc2018_20b2_raster100m/CLC2018_CLC2018_V2018_20b2.tif CLC2018_V2018_20b2.tif 
CLCRaster <- raster(file.path('BASE', 'EGV', 'CLC2018_V2018_20b2.tif'))
## load a reclassificaiton table, reclassify
CLCReclassTable <- read_excel(file.path('data', 'CLC2018_reclassification.xlsx'))
CLCReclassMatrix <- as.matrix(CLCReclassTable[, c('CLASS', 'CLASS.new')])
landcover <- reclassify(CLCRaster, CLCReclassMatrix)
landcover <- ratify(landcover)
# reproject now, landcover should be supposedly lighter
if(projection(landcover) != projection(theCRS)) { # reproject only if needed
  landcover <- projectRaster(from=landcover, to=refRaster, res=theResolution, method='ngb')
} else {
  landcover <- resample(landcover, refRaster, method='ngb')
}
names(landcover) <- "CORINE2018"
rm(CLCRaster)
## calculate percentage maps
landcover.lvl <- data.frame(unique(CLCReclassTable[,c('CLASS.new', 'Description.new')]))
names(landcover.lvl) <- c('ID', 'landcover')
levels(landcover) <- landcover.lvl
landcover <- deratify(landcover, 'landcover')
landcover <- layerize(landcover)
names(landcover) <- landcover.lvl$landcover[1:15] 
# 
filter <- make_circ_filter(theRadius, theResolution)
filterSum <- sum(filter, na.rm=TRUE)
landcoverStack <- stack()
for(n in names(landcover)) {
  pct <- focal(landcover[[n]], filter, fun=function(x) sum(x, na.rm=TRUE)/filterSum)
  landcoverStack <- addLayer(landcoverStack, pct)
}
names(landcoverStack) <- names(landcover)

# ### bonus: add distance from urban areas
# we don't calculate this for Ibex
# urb <- landcover$Urban.areas > 0
# urb[urb==0] <- NA
# urbandistance <- distance(urb)
# names(urbandistance) <- 'urban.distance'

#### add landcover to pack
EGVStack <- addLayer(EGVStack, landcoverStack)

#### cleanup
rm(landcover, landcoverStack, landcover.lvl, pct)

#### add solar radiation -------------------------------------------------------
# using a cropped version of CM-SAF-gh_0_year raster, obtained from https://joint-research-centre.ec.europa.eu/pvgis-photovoltaic-geographical-information-system/pvgis-data-download/cm-saf-solar-radiation_en
RadiationRaw <- raster(file.path('BASE', 'EGV', 'CM-SAF-gh_0_year.tif'))
# reproject and crop
if(projection(RadiationRaw) != projection(theCRS)) { # reproject only if needed
  radiation <- projectRaster(from=RadiationRaw, to=refRaster, res=theResolution, method='bilinear') # bilinear interpolation since it's a continuous variable
} else {
  radiation <- resample(radiation, refRaster, method='ngb')
}
names(radiation) <- "Solar.Radiation"

#### add to pack
EGVStack <- addLayer(EGVStack, radiation)
#### cleanup
rm(RadiationRaw, radiation)

#### add WORLDCLIM data --------------------------------------------------------
# use 30 arc second WORLDCLIM bioclimatic variables
BIORaw <- list.files(path='/srv/archivio/Cartografia/WORLDCLIM_v2/', pattern='wc2.0_bio_30s_..\\.tif', full.names=TRUE)
WORLDCLIMStack <- stack(BIORaw)
names(WORLDCLIMStack) <- paste0("BIO", substr(list.files(path='/srv/archivio/Cartografia/WORLDCLIM_v2/', pattern='wc2.0_bio_30s_..\\.tif'),14,16))
# crop to AOI
AOI_4326 <- st_transform(AOI, st_crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
WCStack <- crop(WORLDCLIMStack, AOI_4326)
# reproject
WCStack <- projectRaster(from=WCStack, to=refRaster, res=theResolution, method='bilinear')

#### add to pack
EGVStack <- addLayer(EGVStack, WCStack)
#### cleanup
rm(WORLDCLIMStack, WCStack, AOI_4326)

## do not use road network, at least ath Alpine Space Area scale, it would be too cumbersome to have all road data.
## do not use distance from inhabited areas, as well.

# shape up final EGV stack, setMinMax takes ages, so do it once and for all
EGVStack[is.na(EGVStack)] <- 0
EGVStack <- setMinMax(EGVStack) # this takes ages

# stop parallelizing
endCluster()

## Not Run:
# to reconstruct from separte TIFFs:
# EGVLayers <- list.files(path='BASE/EGV/', pattern='EGV_layer.*\\.tif', full.names=TRUE)
# EGVStack <- stack(EGVLayers)
# names(EGVStack) <- c("Elevation",  "Aspect.flat", "Aspect.N", "Aspect.NE", "Aspect.E", "Aspect.SE", "Aspect.S", "Aspect.SW", "Aspect.W", "Aspect.NW",  "Slope.level", "Slope.moderate", "Slope.strong", "Slope.steep", "Roughness", "Urban.areas", "Agriculture", "Seminatural.vegetation", "Woodlands", "Natural.grasslands", "Moors.and.heathland", "Shrubland", "Beaches.dunes.sands", "Bare.rocks", "Sparsely.vegetated.areas", "Burnt.areas", "Glaciers.and.perpetual.snow", "Marshes.and.bogs", "Salt.marshes", "Water.bodies", "Solar.Radiation", "BIO_01", "BIO_02", "BIO_03", "BIO_04", "BIO_05", "BIO_06", "BIO_07", "BIO_08", "BIO_0data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==9", "BIO_10", "BIO_11", "BIO_12", "BIO_13", "BIO_14", "BIO_15", "BIO_16", "BIO_17", "BIO_18", "BIO_19")  

## Not Run:
# to check whether everything is OK (NAs, values etc.)
# for(s in names(EGVStack)) {cat(s,'\n'); print(EGVStack[[s]])}

#### save to disk --------------------------------------------------------------
writeRaster(EGVStack, filename=paste0('./BASE/EGV/', 'EGVStack.grd'), overwrite=TRUE, format='raster')
# read in with  EGVStack <- raster(filename=paste0('./BASE/EGV/', EGVStackName))

# save as ERDAS IMAGINE
writeRaster(EGVStack, filename=paste0('./BASE/EGV/', "EGV.img"), overwrite=TRUE, format='HFA', options="COMPRESSED=YES")

# also save in a more "neutral" format, preserving 'band names'...
for(n in names(EGVStack)) {
  writeRaster(EGVStack[[n]], filename=paste0('./BASE/EGV/EGV_', n, '.tif'), format='GTiff', overwrite=TRUE)
}

#### extract a subset for Lombardy Alpine space --------------------------------
AOI.lombardy <- st_read(dsn='BASE', layer='AOI Lombardia 3035') # already in EPSG:3035, no need reprojecting
EGVStack.lombardy <- terra::crop(EGVStack, AOI.lombardy, snap="out")
writeRaster(EGVStack.lombardy, filename=paste0('./BASE/EGV/', 'EGV-Lombardy.grd'), overwrite=TRUE, format='raster')
writeRaster(EGVStack.lombardy, filename=paste0('./BASE/EGV/', "EGV-Lombardy.img"), overwrite=TRUE, format='HFA', options="COMPRESSED=YES")


### End Of File ###
