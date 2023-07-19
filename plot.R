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

library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(rnaturalearth)
library(sp)
library(ggspatial)
library(vegan)

# ----------------------------------------------------------


# ------------- CREATE PLOT DIRECTORY ----------------------
# create 'plots' directory to put the plots in

createDirectory <- function(){

	if(!file.exists('Output/plots')){
		dir.create(path='Output/plots')
	}

}

# ----------------------------------------------------------


# ------------------- READ REPLICATES ----------------------
# read all the replicate runs in the output directory

readReplicates <- function(){

	# get list of replicate files
	files <- list.files(path='Output',pattern='*.csv',full.names=TRUE)
	# create container for simulation output
	replicateContainer <- matrix(ncol=4,nrow=0) %>% as.data.frame()
	names(replicateContainer) <- c('lon','lat','diversity','rep')

	# get diversity matrix output switch from simulation specifics file
	PAR_file <- file('Output/PARAM_file.txt',open='r')
	on.exit(close(PAR_file))
	PAR_file_lines <- readLines(PAR_file)
	spec_line <- PAR_file_lines[which(startsWith(PAR_file_lines,'Print complete diversity matrix to output'))]
	spec <- str_split(spec_line,pattern=': ')[[1]][2]
	if(spec == 'active'){
		spec_switch <- TRUE
	} else {
		spec_switch <- FALSE
	}

	for(i in 1:length(files)){
		dat <- read.csv(files[i],header=TRUE,row.names=1)
		# add replicate to the container
		if(spec_switch){
			# calculate diversity
			dat <- dat %>% dplyr::mutate(diversity=dat[,-c(1,2)] %>% vegan::decostand(method='pa') %>% rowSums()) %>%
							dplyr::select(c(lon,lat,diversity)) %>% dplyr::mutate(rep=i)

			# add to container
			replicateContainer <- rbind(replicateContainer,dat)

		} else {
			replicateContainer <- rbind(replicateContainer,dat %>% dplyr::mutate(rep=i))
		}
	}

	# get the summary of replicates
	df_sum <- replicateContainer %>% dplyr::group_by(lat,lon) %>% dplyr::summarise(simulation_mean=mean(diversity),simulation_sd=sd(diversity)) %>% 
										dplyr::ungroup() %>% dplyr::mutate(simulation_mean=simulation_mean/max(simulation_mean,na.rm=TRUE))
	# transform zeros to NA (inactive cells)
	df_sum$simulation_mean[df_sum$simulation_mean==0] <- NA
	df_sum$simulation_sd[df_sum$simulation_sd==0] <- NA

	return(df_sum)
	
}

# ----------------------------------------------------------


# --------------- COMBINE SIMULATION WITH DATA -------------

createPlotDF <- function(){

	# get simulation grid file from simulation specifics file
	PAR_file <- file('Output/PARAM_file.txt',open='r')
	on.exit(close(PAR_file))
	PAR_file_lines <- readLines(PAR_file)
	grid_line <- PAR_file_lines[which(startsWith(PAR_file_lines,'\t- grid:'))]
	grid_file <- str_split(grid_line,pattern=': ')[[1]][2]

	dat_emp <- read.csv(paste0('GridFiles/',grid_file),header=TRUE,row.names=1) %>% 
				dplyr::select(x,y,diversity)

	# get simulation run
	dat_sim <- readReplicates()

	# combine empirical data with simulation output
	dat_sim_emp <- dat_sim %>% dplyr::inner_join(.,dat_emp,by=c('lon'='x','lat'='y')) %>%
					dplyr::mutate(residual=(simulation_mean - diversity))

	return(dat_sim_emp)	

}

# ----------------------------------------------------------


# --------------- CREATE BACKGROUND SHAPEFILES -------------

getSHP <- function(switch){

	# Contiguous USA for NABB
	if(switch=='NABB'){

		# get contiguous USA shapefile
		USA <- rnaturalearth::ne_countries(country='united states of america',scale='medium',returnclass='sf') %>% sf::st_geometry() %>%
						sf::st_crop(xmin=-130,xmax=-50,ymin=20,ymax=55)
		# transform to Albers projection (equal area)
		shp <- USA %>% sf::st_transform(crs=5070)

	} 

	# North Atlantic for CPR
	else if(switch=='CPR'){

		# define the Lambert Conformal Conic projection
		crs_atl <- CRS('+proj=lcc +lon_0=-20 +lat_1=40 +lat_2=60')

		ex <- sf::st_sf(geom=sf::st_sfc(sf::st_point(c(-2200000,5436156)),sf::st_point(c(-2200000,7500000)),
										sf::st_point(c(1850000,5436156)),sf::st_point(c(1850000,7500000)))) %>% 
								sf::st_bbox() %>% raster::extent()

		# create the simulation grid raster
		r <- raster::raster(ncol=45,nrow=14,crs=crs_atl,ext=ex)

		# create land layer for plotting
		shp_land <- rnaturalearth::ne_countries(scale='medium') %>% sf::st_as_sf() %>% sf::st_geometry() %>%
						sf::st_crop(., sf::st_sf(geom=sf::st_sfc(sf::st_point(c(-65,33)),sf::st_point(c(-65,70)),
											sf::st_point(c(20,33)),sf::st_point(c(20,70)))) %>% sf::st_bbox()) %>%
						sf::st_transform(crs=crs_atl)
		# crop land layer to raster size
		shp <- shp_land %>% sf::st_crop(r) 

	} else {

		print('switch unrecognised')
		shp <- NULL

	}

	return(shp)

}

# ----------------------------------------------------------


# ---------------------- MAKE PLOTS ------------------------

plotSimulation <- function(){

	createDirectory()

	# get the correct switch value
	PAR_file <- file('Output/PARAM_file.txt',open='r')
	on.exit(close(PAR_file))
	PAR_file_lines <- readLines(PAR_file)
	grid_line <- PAR_file_lines[which(startsWith(PAR_file_lines,'\t- grid:'))]
	grid_file <- str_split(grid_line,pattern=': ')[[1]][2]
	switch <- str_split(grid_file,pattern='_')[[1]][1]

	# get the correct shapefile for background
	shp <- getSHP(switch)

	# get the plot dataframe
	df_plot <- createPlotDF()

	# get minimal value
	MinVal <- min(c(df_plot$diversity,df_plot$simulation_mean),na.rm=TRUE)

	# plot the simulation mean
	p <- ggplot() + geom_sf(data=shp) + geom_tile(data=df_plot,aes(x=lon,y=lat,fill=simulation_mean),color='black') + theme_bw() +
			scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='simulation mean',limits=c(MinVal,1)) +
			ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
				style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
			ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
			theme(axis.title=element_blank(), legend.title=element_text(size=10))

	p %>% ggsave('Output/plots/simulation_mean.png',.,device='png',width=20,height=12,units='cm')

	# plot the simulation standard deviation
	p <- ggplot() + geom_sf(data=shp) + geom_tile(data=df_plot,aes(x=lon,y=lat,fill=simulation_sd),color='black') + theme_bw() +
			scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='simulation sd') +
			ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
				style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
			ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
			theme(axis.title=element_blank(), legend.title=element_text(size=10))

	p %>% ggsave('Output/plots/simulation_sd.png',.,device='png',width=20,height=12,units='cm')

	# plot the simulation residuals map
	p <- ggplot() + geom_sf(data=shp) + geom_tile(data=df_plot,aes(x=lon,y=lat,fill=residual),color='black') + theme_bw() +
			scale_fill_gradient2(low=rgb(1,0,0,0.65),high=rgb(0,1,0,0.65),mid=rgb(1,1,1,0.65),
									name='residuals',limits=c(-1,1),na.value=rgb(0.5,0.5,0.5,0.65)) +
			ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
				style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
			ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
			theme(axis.title=element_blank(), legend.title=element_text(size=10))

	p %>% ggsave('Output/plots/residual_map.png',.,device='png',width=20,height=12,units='cm')

	# plot the simulation residual distribution
	p <- df_plot %>% ggplot() + geom_histogram(aes(residual)) + theme_bw() 
	p %>% ggsave('Output/plots/residual_distribution.png',.,device='png',width=10,height=10,units='cm')

	# calculate the regression coefficients between simulation and empirical data
	model <- lm(data=df_plot,simulation_mean~diversity)
	coeff <- summary(model)$coefficients %>% as.data.frame() %>% dplyr::select(Estimate)
	# calculate the Pearson correlation coefficient
	correl <- cor.test(df_plot$diversity,df_plot$simulation_mean,method='pearson')$estimate %>% as.numeric()

	# plot the scatter
	p <- df_plot %>% ggplot() + geom_point(aes(x=diversity,y=simulation_mean),pch=19) + 
						geom_abline(intercept=coeff[1,1], slope=coeff[2,1], color='red', linewidth=1) +
						geom_abline(intercept=0,slope=1,linetype=2,linewidth=0.75) +
						xlab('data') + ylab('simulation') +
						theme_bw() + xlim(0,1) + ylim(0,1) +
						labs(title=paste0('r = ',round(correl,4)))

	p %>% ggsave('Output/plots/scatter.png',.,device='png',width=10,height=10,units='cm')

}

# ----------------------------------------------------------

# make the plots
plotSimulation()