library(sp)
library(rgdal)
library(raster)
library(tidyverse)
library(vegan)

if(F){
dat <- read.csv('Output/grid_0000.csv',header=T,row.names=1)

datpa <- decostand(dat[,3:ncol(dat)],method='pa')
count <- rowSums(datpa) 

r <- raster(ncol=45,nrow=14,xmn=-180,xmx=180,ymn=-90,ymx=90)
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

test <- as.data.frame(rasterToPoints(r))
test$z <- as.vector(count)

coordinates(test) <- ~x+y
proj4string(test)=CRS("+init=epsg:4326")

gridded(test) <- TRUE

r <- raster(test)
plot(r)
}

step_list <- rep(0,1001)
total_list <- rep(0,1001)

for(i in 0:1000){

	dat <- read.csv(paste0('Output/grid_',sprintf("%04d", i),'.csv'),header=T,row.names=1)
	step_list[i+1] <- i*50000
	total_list[i+1] <- ncol(dat)-2
}

df_plot <- data.frame(time=step_list, div=total_list)

p <- df_plot %>% ggplot() + geom_line(aes(x=time,y=div),linewidth=1) + theme_bw() +
					xlab('turnover') + ylab('global diversity')

p %>% ggsave('global_div.png',.,device='png',width=15,height=10,units='cm')