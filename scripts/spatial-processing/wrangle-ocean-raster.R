# turn ocean netcdf into tif and save

library(terra)

dat <- rast('data/ocean/AQUA_MODIS.20230301_20230331.L3m.MO.SST.sst.9km.nc')
plot(dat)

writeRaster(dat[[1]], 'data/ocean/sst_raster.tif')
