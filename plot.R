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

}

if(F){

dat <- read.csv('Output/grid_0000.csv',header=T,row.names=1)
datLatLon <- dat %>% select(c(1,2))
datPA <- dat %>% select(-c(1,2)) %>% decostand(.,method='pa') %>% mutate(total=rowSums(.)) %>% pull(total)
datLatLon <- datLatLon %>% mutate(total=datPA) %>% group_by(lat) %>% summarise(aver=mean(total))

print(datLatLon)

g <- datLatLon %>% ggplot() + geom_line(aes(x=lat,y=aver),linewidth=1) +  coord_flip() +
					theme_bw() + labs(x = 'latitude (°)', y = 'mean local diversity') +
					ylim(0,7) + scale_x_continuous(breaks = c(-90,-45,0,45,90), limits=c(-95,95))

g %>% ggsave('timeseries.png',.,device='png',height=15,width=10,units='cm')

}

if(F){

r <- raster(ncol=45,nrow=14,xmn=-180,xmx=180,ymn=-90,ymx=90)
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

test <- as.data.frame(rasterToPoints(r))

dat <- read.csv('grid_0000.csv',header=T,row.names=1)

datpa <- decostand(dat[,3:ncol(dat)],method='pa')
count <- rowSums(datpa) 

test <- as.data.frame(rasterToPoints(r))
test$z <- as.vector(count)
coordinates(test) <- ~x+y
proj4string(test)=CRS("+init=epsg:4326")

gridded(test) <- TRUE

r2 <- raster(test)

test_spdf <- as(r2, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

g <- test_df %>% ggplot() + geom_tile(aes(x=x,y=y,fill=value)) + theme_bw() + 
		scale_x_continuous(breaks=seq(-180,180,45), limits=c(-185,180)) + scale_y_continuous(breaks=seq(-90,90,45),limits=c(-95,95)) +
		labs(x='longitude',y='latitude') + 
		scale_fill_viridis_c(name='local diversity',option='magma',breaks=seq(1,16,3)) + 
		theme(panel.ontop=T,panel.background = element_rect(fill = NA))

g %>% ggsave('raster.png',.,device='png',width=15,height=10,units='cm')

}

if(T){

dat <- read.csv('Output/grid_0000.csv',header=T,row.names=1)

dat <- dat %>% select(c(1,2,2711,floor((ncol(.)-2)/3),floor((ncol(.)-2)*2/3),ncol(.))) 

r <- raster(ncol=45,nrow=14,xmn=-180,xmx=180,ymn=-90,ymx=90)
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

test_df <- data.frame(matrix(ncol=4,nrow=0))
colnames(test_df) <- c('count','x','y','species')

test <- as.data.frame(rasterToPoints(r))

for(i in 3:ncol(dat)){

	test <- as.data.frame(rasterToPoints(r))
	test$count <- as.vector(dat[,i])
	coordinates(test) <- ~x+y
	proj4string(test)=CRS("+init=epsg:4326")

	gridded(test) <- TRUE

	r2 <- raster(test)

	test_spdf <- as(r2, "SpatialPixelsDataFrame")
	test_df_sub <- as.data.frame(test_spdf) %>% mutate(species=colnames(dat)[i])
	colnames(test_df_sub) <- c('count','x','y','species')

	test_df <- rbind(test_df,test_df_sub)

}

test_df$count[test_df$count==0] <- NA

p <- test_df %>% ggplot() + geom_tile(aes(x=x,y=y,fill=count),col='black',linewidth=0.5) + facet_wrap(~species) + theme_bw() + 
		scale_x_continuous(breaks=seq(-180,180,45), limits=c(-185,180)) + scale_y_continuous(breaks=seq(-90,90,45),limits=c(-95,95)) +
		theme(panel.ontop=T,panel.background = element_rect(fill = NA),
				panel.grid.major = element_line(color = "grey",linewidth = 0.75,linetype = 2)) + 
		scale_fill_viridis_c(name='population count',na.value='white',breaks=seq(1,16,3), option='magma') + 
		labs(x='longitude',y='latitude') 

p %>% ggsave('grid_species.png',.,width=30,height=20,units='cm',device='png')

}

if(F){

dat <- read.csv('grid_0000.csv',header=T,row.names=1) 

dat_spec <- dat %>% select(-c(1,2))

dat_spec.h <- dat_spec %>% decostand(.,method='hellinger')

PCA <- rda(dat_spec.h)
Evalues <- eigenvals(PCA) 
lambda <- sprintf("%.2f",round(Evalues[1:2] / sum(Evalues) *100, 2))

site.scores <- scores(PCA, scaling=1, display="sites")
df_plot <- site.scores %>% as.data.frame() %>% mutate(lat = dat %>% pull(lat) %>% abs())

p <- df_plot %>% ggplot() + geom_point(aes(x=PC1,y=PC2,col=lat),size=1) + xlab(paste0('PC1 (',lambda[1],' %)')) + ylab(paste0('PC2 (',lambda[2],' %)')) + 
				theme_bw() + scale_colour_viridis_c(option='magma',name='latitude') + 
				geom_hline(yintercept=0,linetype=2) + geom_vline(xintercept=0,linetype=2)

p %>% ggsave('PCA_output.png',.,device='png',width=15,height=10,units='cm')

}