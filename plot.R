library(sp)
library(rgdal)
library(raster)

dat <- read.csv('grid.csv',header=T,row.names=1)

r <- raster(ncol=45,nrow=14,xmn=-180,xmx=180,ymn=-90,ymx=90)
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

test <- as.data.frame(rasterToPoints(r))
test$z <- dat[,3]

coordinates(test) <- ~x+y
proj4string(test)=CRS("+init=epsg:4326")

gridded(test) <- TRUE

r <- raster(test)
plot(r)