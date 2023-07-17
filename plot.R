# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Simulation of Mechanistic Model
#####################################
# PLOT simulation
#####################################


# --------------------- LOAD PACKAGES ----------------------

library(sp)
library(rgdal)
library(raster)
library(tidyverse)
library(vegan)
library(gganimate)
library(gifski)
library(transformr)

# ----------------------------------------------------------


# ---------------------- CREATE GRID -----------------------
# the grid from Worm and Tittensor (2018) is not a connected grid in WGS84
# for plotting the results, create a new global grid.

createGrid <- function(){

	r <- raster::raster(ncol=45,nrow=14,xmn=-180,xmx=180,ymn=-90,ymx=90)

	df <- r %>% as.data.frame(xy=TRUE)

	return(df %>% dplyr::select(x,y))

}

# ----------------------------------------------------------


# ------------- COLLECT TEMPORAL INFORMATION ---------------
# extract information on the time evolution of the simulation

temporalSimulation <- function(){

	# list all csv output files
	list_files <- list.files(path='Output',pattern='*.csv',full.names=TRUE)

	# create container for longitude / latitude diversity
	datLonLat <- matrix(nrow=0,ncol=4) %>% as.data.frame()
	names(datLonLat) <- c('lon','lat','div','time')
	# create container for total diversity
	datTot <- matrix(nrow=0,ncol=2) %>% as.data.frame()
	names(datTot) <- c('time','div')

	# get print interval from simulation specifics file
	PAR_file <- file('Output/PARAM_file.txt',open='r')
	on.exit(close(PAR_file))
	PAR_file_lines <- readLines(PAR_file)
	turn_line <- PAR_file_lines[which(startsWith(PAR_file_lines,'Number of turnovers after which'))]
	Nturn <- str_split(turn_line,pattern=': ')[[1]][2] %>% as.numeric()

	k <- 0
	# for each printed grid
	for(file in list_files){
		cat('\r','file ',k,' of ',length(list_files))
		dat <- read.csv(file,header=TRUE,row.names=1)
		# get the number of print
		i <- str_split(file,pattern='_')[[1]][2] %>% str_split(pattern='[.]')
		i <- i[[1]][1] %>% as.numeric()

		# calculate diversity for each grid cell
		dat <- dat %>% dplyr::mutate(div=dat[,-c(1,2)] %>% vegan::decostand(method='pa') %>% rowSums())

		# add to the container and add the simulation time = number of turnovers
		datLonLat <- rbind(datLonLat,
								dat %>% dplyr::select(c(lon,lat,div)) %>% dplyr::mutate(time=i*Nturn))

		# add total diversity to the container
		datTot <- rbind(datTot,data.frame(time=i*Nturn,div=ncol(dat)-2))
		k <- k + 1

	}

	return(list(datLonLat,datTot))

}

# ----------------------------------------------------------


# ------------- CREATE PLOT DIRECTORY ----------------------
# create 'plots' directory to put the plots in

createDirectory <- function(){

	if(!file.exists('Output/plots')){
		dir.create(path='Output/plots')
	}

}


# --------- PLOT GLOBAL DIVERSITY EVOLUTION ----------------

globalDiversity <- function(datTot){

	createDirectory()

	p <- datTot %>% ggplot() + geom_line(aes(x=time,y=div)) + theme_bw() +
						xlab('turnover') + ylab('global number of species')

	p %>% ggsave('Output/plots/global_diversity.png',.,width=15,height=10,units='cm')

}

# ----------------------------------------------------------


# ---------- PLOT LOCAL DIVERSITY EVOLUTION ----------------

localDiversityGrid <- function(coord,datLonLat){

	createDirectory()

	# replace coordinates with newly created coordinates for correct plotting
	datLonLat$lon <- coord$x
	datLonLat$lat <- coord$y

	g <- datLonLat %>% ggplot() + geom_tile(aes(x=lon,y=lat,fill=div)) + theme_bw() + 
			#scale_x_continuous(breaks=seq(-180,180,45), limits=c(-185,180)) + scale_y_continuous(breaks=seq(-90,90,45),limits=c(-95,95)) +
			labs(title='turnover: {as.integer(frame_time)}',x='longitude',y='latitude') + 
			scale_fill_viridis_c(name='local diversity',option='magma') + 
			theme(panel.ontop=T,panel.background = element_rect(fill = NA)) + transition_time(time) + ease_aes('linear')

	# might need some adjustment depending on the simulation:
	# res = resolution
	# fps = frames per second
	g %>% gganimate::animate(width=15,height=10,units='cm',nframes=length(unique(datLonLat$time)),res=150,fps=10,renderer = gifski_renderer()) %>% 
			gganimate::anim_save('Output/plots/timeseries_grid.gif',.)

}

# ----------------------------------------------------------


# ----- PLOT LATITUDINAL LOCAL DIVERSITY EVOLUTION ---------

localDiversityLat <- function(datLonLat){

	createDirectory()

	datLat <- datLonLat %>% dplyr::group_by(lat,time) %>% dplyr::summarise(div=mean(div))

	g <- datLat %>% ggplot() + geom_line(aes(x=lat,y=div),linewidth=1) +  coord_flip() +
					theme_bw() + labs(title = 'turnover: {as.integer(frame_time)}', x = 'latitude (°)', y = 'mean local diversity') +
					ylim(0,max(datLat$div)+1) + scale_x_continuous(breaks = c(-90,-45,0,45,90), limits=c(-95,95)) +
					transition_time(time) + ease_aes('linear') 

	g %>% gganimate::animate(width = 10, height=15,units='cm',res=150,nframes=length(unique(datLonLat$time)),fps=10,renderer = gifski_renderer()) %>% 
				gganimate::anim_save('Output/plots/timeseries_lat.gif',.)

}

# ----------------------------------------------------------


# --------------- PLOT SPECIES GRID ------------------------
# plot number of populations of species 'index' in grid
# the last printed output is taken
# if index = -1; the species that entered last is taken
# if index = 1; the species that entered first is taken

speciesGrid <- function(coord,index){

	createDirectory()

	# read the final output file
	list_files <- list.files(path='Output',pattern='*.csv',full.names=TRUE)
	dat <- read.csv(list_files[length(list_files)],header=TRUE,row.names=1)
	# replace zeros by NA
	dat[dat==0] <- NA
	# replace coordinates for plotting
	dat$lon <- coord$x
	dat$lat <- coord$y

	# get the correct species
	if(index == -1){
		i <- ncol(dat)
	}
	if(index == 1){
		i <- 3
	}

	# get number of local populations from PARAM file
	PAR_file <- file('Output/PARAM_file.txt',open='r')
	on.exit(close(PAR_file))
	PAR_file_lines <- readLines(PAR_file)
	loc_line <- PAR_file_lines[which(startsWith(PAR_file_lines,'Number of local populations'))]
	Nloc <- str_split(loc_line,pattern=': ')[[1]][2] %>% as.numeric()

	p <- dat %>% ggplot() + geom_tile(aes(x=lon,y=lat,fill=dat[,i]),col='black',linewidth=0.5)  + theme_bw() + 
		#scale_x_continuous(breaks=seq(-180,180,45), limits=c(-185,180)) + scale_y_continuous(breaks=seq(-90,90,45),limits=c(-95,95)) +
		theme(panel.ontop=T,panel.background = element_rect(fill = NA)) + 
		scale_fill_viridis_c(name='population count',na.value='white',breaks=seq(1,Nloc,3), option='magma') + 
		labs(x='longitude',y='latitude') 

	p %>% ggsave(paste0('Output/plots/species_',index,'_grid.png'),.,width=18,height=10,units='cm',device='png')

}

# ----------------------------------------------------------


# --------------------- PLOT PCA ---------------------------
# plot a PCA of the final output file

plotPCA <- function(){

	createDirectory()

	# get the final output file
	list_files <- list.files(path='Output',pattern='*.csv',full.names=TRUE)
	dat <- read.csv(list_files[length(list_files)],header=TRUE,row.names=1) 

	# get diversity matrix
	dat_spec <- dat %>% select(-c(1,2))
	# Hellinger transformation
	dat_spec.h <- dat_spec %>% vegan::decostand(.,method='hellinger')
	# calculate the PCA
	PCA <- vegan::rda(dat_spec.h)
	# get cell scores and add latitude of cells
	site.scores <- vegan::scores(PCA, scaling=1, display="sites")
	df_plot <- site.scores %>% as.data.frame() %>% dplyr::mutate(lat = dat %>% dplyr::pull(lat) %>% abs())
	# get proportion explained
	percs <- summary(PCA)$cont$importance[2,1:2]*100 
	percs <- percs %>% round(.,digits=2)

	p <- df_plot %>% ggplot() + geom_point(aes(x=PC1,y=PC2,col=lat),size=1) + 
					xlab(paste0('PC1 (',percs[1],' %)')) + ylab(paste0('PC2 (',percs[2],' %)')) + 
					theme_bw() + scale_colour_viridis_c(option='magma',name='latitude') + 
					geom_hline(yintercept=0,linetype=2) + geom_vline(xintercept=0,linetype=2)

	p %>% ggsave('Output/plots/PCA_output.png',.,device='png',width=15,height=10,units='cm')

}

# ----------------------------------------------------------


# ---------------- RUN ALL THE PLOTS -----------------------
t1 <- Sys.time()

coord <- createGrid()
res <- temporalSimulation()
datLonLat <- res[[1]]
datTot <- res[[2]]

globalDiversity(datTot)
localDiversityGrid(coord,datLonLat)
localDiversityLat(datLonLat)
speciesGrid(coord,-1)
speciesGrid(coord,1)
plotPCA()

t2 <- Sys.time()

td <- seconds_to_period(difftime(t2,t1,units='secs'))
cat('\n')
print(paste0('wall time: ',sprintf('%02d:%02d:%02d',td@hour, minute(td), round(second(td))) ))