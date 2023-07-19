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
library(vegan)
library(ggpubr)

# ----------------------------------------------------------


# ------------- CREATE PLOT DIRECTORY ----------------------
# create 'plots' directory to put the plots in

createDirectory <- function(){

	if(!file.exists('Output/plots')){
		dir.create(path='Output/plots')
	}

}

# ----------------------------------------------------------


# --------------------- SWITCH CHECK -----------------------

readPARAMs <- function(){

	# get simulation specifics file
	PAR_file <- file('Output/PARAM_file.txt',open='r')
	on.exit(close(PAR_file))
	PAR_file_lines <- readLines(PAR_file)
	# get the simulation grid file
	grid_line <- PAR_file_lines[which(startsWith(PAR_file_lines,'\t- grid:'))]
	grid_file <- str_split(grid_line,pattern=': ')[[1]][2]
	# get the habitat area switch
	hab_area <- PAR_file_lines[which(startsWith(PAR_file_lines,'Habitat area scaling'))]
	hab_area <- str_split(hab_area,pattern=': ')[[1]][2]
	if(hab_area == 'active'){
		hab_switch <- TRUE
	} else {
		hab_switch <- FALSE
	}
	# get the temperature switch
	temp <- PAR_file_lines[which(startsWith(PAR_file_lines,'Temperature-dependent speciation rates'))]
	temp <- str_split(temp,pattern=': ')[[1]][2]
	if(temp == 'active'){
		temp_switch <- TRUE
	} else {
		temp_switch <- FALSE
	}
	spec <- PAR_file_lines[which(startsWith(PAR_file_lines,'Print complete diversity matrix to output'))]
	spec <- str_split(spec,pattern=': ')[[1]][2]
	if(spec == 'active'){
		spec_switch <- TRUE
	} else {
		spec_switch <- FALSE
	}
	# get number of local populations
	Nloc <- PAR_file_lines[which(startsWith(PAR_file_lines,'Base number of local populations'))]
	Nloc <- str_split(Nloc,pattern=': ')[[1]][2] %>% as.numeric()
	# get longitudinal dimension
	Nlon <- PAR_file_lines[which(startsWith(PAR_file_lines,'Number of longitudinal cells'))]
	Nlon <- str_split(Nlon,pattern=': ')[[1]][2] %>% as.numeric()
	# get latitudinal dimension
	Nlat <- PAR_file_lines[which(startsWith(PAR_file_lines,'Number of latitudinal cells'))]
	Nlat <- str_split(Nlat,pattern=': ')[[1]][2] %>% as.numeric()

	return(list(grid_file,hab_switch,temp_switch,spec_switch,Nloc,Nlon,Nlat))

}

# ----------------------------------------------------------


# ------------ CREATE DIVERSITY MATRIX ---------------------
# create the diversity matrix from the list of species identifiers in the 'Output' directory

getSpeciesTable <- function(){

	# get list of replicate run Output files
	files <- list.files(path='Output', pattern="grid_*", full.names=TRUE, recursive=FALSE)

	# get simulation parameters
	PARs <- readPARAMs()
	grid_file <- PARs[[1]]
	hab_switch <- PARs[[2]]
	Nloc <- PARs[[5]]
	
	# for each replicate run
	for(file in files){

		# get the list of species identifiers
		ind_list <- read.csv(file,header=TRUE,row.names=1)
		# get the simulation grid
		grid_list <- read.csv(paste0('GridFiles/',grid_file),header=TRUE,row.names=1) 

		# calculate the local number of populations
		if(hab_switch){
			grid_list <- grid_list %>% dplyr::mutate(Nloc_i=floor(Nloc*grad))
		} else {
			grid_list <- grid_list %>% dplyr::mutate(Nloc_i=Nloc)
		}

		# get list of species in the simulation
		spec_list <- sort(unique(ind_list$species))

		# create a container for final diversity grid
		df_out <- matrix(nrow=nrow(grid_list),ncol=2+length(spec_list)) %>% as.data.frame()
		names(df_out) <- c('lon','lat',spec_list)

		# counter in ind_list
		count <- 1
		# for each grid cell
		for(i in 1:nrow(grid_list)){
			# get the populations that are present in the cell
			ind_sub <- ind_list[count:(count+grid_list$Nloc_i[i]-1),] %>% as.data.frame() %>% dplyr::rename(species='.')
			# calculate the abundance of each species
			ind_sub_sum <- ind_sub %>% dplyr::group_by(species) %>% dplyr::summarise(count=length(species))
			# create a container for the species abundances in the cell
			spec_container <- data.frame(species=spec_list,count=0)
			# order the abundances, so non-present species have 0 abundance
			for(j in 1:nrow(ind_sub_sum)){
				spec_container$count[spec_container$species==ind_sub_sum$species[j]] <- ind_sub_sum$count[j]
			}
			# add the cell to the final container
			df_out[i,] <- c(grid_list$lon[i],grid_list$lat[i],spec_container$count) %>% t()

			# the counter increases with the number of populations in the cell
			count <- count + grid_list$Nloc_i[i]
		}

		# get the file identifier
		file_number <- str_split(file,pattern='_')[[1]][2]
		file_number <- str_split(file_number,pattern='[.]')[[1]][1]

		#  write out the diversity matrix for this replicate run
		write.csv(df_out,paste0('Output/spec_',file_number,'.csv'))
	}

}

# ----------------------------------------------------------


# ------------ COMBINE SIMULATION WITH DATA ----------------
# make a dataframe combining the simulation with the data

readReplicatesDiversity <- function(){

	# get simulation parameters
	PARs <- readPARAMs()
	grid_file <- PARs[[1]]
	spec_switch <- PARs[[4]]

	# calculate the diversity matrices
	if(spec_switch){
		getSpeciesTable()
	}

	# get list of replicate run diversity matrix files
	if(spec_switch){
		files <- list.files(path='Output', pattern="spec_*", full.names=TRUE, recursive=FALSE)
	} else {
		files <- list.files(path='Output', pattern="grid_*", full.names=TRUE, recursive=FALSE)
	}

	# get the empirical simulation grid
	dat_emp <- read.csv(paste0('GridFiles/',grid_file),header=TRUE,row.names=1)

	# create container to combine all replicate runs
	replicateContainer <- matrix(ncol=4,nrow=0) %>% as.data.frame()
	names(replicateContainer) <- c('lon','lat','div','rep')

	# combine replicate runs
	for(i in 1:length(files)){

		dat_sim <- read.csv(files[i],header=TRUE,row.names=1) 
		if(spec_switch){
			dat_sim <- dat_sim %>% dplyr::mutate(div=dat_sim[,-c(1,2)] %>% vegan::decostand(method='pa') %>% rowSums()) %>%
						dplyr::select(c(lon,lat,div)) %>% dplyr::mutate(rep=i)
		} else {
			dat_sim <- dat_sim %>% dplyr::mutate(rep=i) %>% dplyr::rename(div=diversity)
		}

		replicateContainer <- rbind(replicateContainer,dat_sim)
	}

	# calculate simulation statistics
	df_sum <- replicateContainer %>% dplyr::group_by(lat,lon) %>% dplyr::summarise(simulation_mean=mean(div),simulation_sd=sd(div)) %>% 
										dplyr::ungroup() %>% 
										dplyr::inner_join(.,dat_emp,by=c('lon'='lon','lat'='lat')) %>% 
										dplyr::mutate(boltz=exp(-0.65 / 8.617e-5 / temp),gradNorm=grad/max(grad)) %>%
										dplyr::mutate(boltzNorm=boltz/max(boltz)) %>%
										dplyr::mutate(simulation_mean=simulation_mean/max(simulation_mean,na.rm=TRUE),
														residualGrad=(simulation_mean - gradNorm),
														residualTemp=(simulation_mean - boltzNorm))

	return(df_sum)
	
}

# ----------------------------------------------------------


# --------------- CALCULATE DISTANCE -----------------------
# calculate the distance between two coordinates in the simulation grid
# take into account periodic boundary conditions

distancePeriodic <- function(x,y){

	# get simulation parameters
	PARs <- readPARAMs()
	Nlon <- PARs[[6]]
	Nlat <- PARs[[7]]

	# calculate distances
	dist_lon <- min(abs(x$lon-y$lon),Nlon - abs(x$lon-y$lon))
	dist_lat <- min(abs(x$lat-y$lat),Nlat - abs(x$lat-y$lat))
	
	dist_tot <- sqrt(dist_lon*dist_lon + dist_lat*dist_lat)

	return(dist_tot)

}

# ----------------------------------------------------------


# ------------------ CALCULATE PCF -------------------------
# calculate the Pair Correlation Function for all replicates

calcPCF <- function(){
	
	# get list of all replicate diversity matrices
	files <- list.files(path='Output', pattern="spec_*", full.names=TRUE, recursive=FALSE)
	
	# create container for PCF calculations
	df_PCF <- matrix(nrow=0,ncol=2) %>% as.data.frame()
	names(df_PCF) <- c('distance','PCF')

	print('PCF calculations')
	# for each replicate
	for(i in 1:length(files)){
		print(paste0('replicate ',i,' of ',length(files)))
		# get the diversity matrix
		dat <- read.csv(files[i],header=TRUE,row.names=1) 

		# for each cell
		for(j in 1:(nrow(dat)-1)){
			# for each other cell
			for(k in (j+1):nrow(dat)){
				# get the pair of cells
				dat_sub <-  dat[c(j,k),]
				# calculate the distance between the pairs
				dist_sub <- dat_sub %>% dplyr::select(c(lon,lat))
				dist <- distancePeriodic(dist_sub[1,],dist_sub[2,])
				if(dist==0){
					print(paste(j,k))
					stop()
				}

				# calculate the empirical PCF
				dat_spec <- dat_sub[,-c(1:2)] %>% t() %>% as.data.frame()
				names(dat_spec) <- c('A','B')
				dat_spec$prod <- dat_spec$A * dat_spec$B
				covar <- mean(dat_spec$prod)
				sumA <- mean(dat_spec$A)
				sumB <- mean(dat_spec$B)
				pcf <- covar/sumA/sumB

				# add the distance and PCF to the container
				df_PCF <- rbind(df_PCF,data.frame(distance=dist,PCF=pcf))
			}
		}

	}
	
	return(df_PCF)

}

# ----------------------------------------------------------


# ------------------ THEORETICAL PCF -----------------------

PCFparam <- function(x,par){
	return(1+1/2/pi*(par[1]/par[2])^2*besselK(x/par[2],0))
}

# ----------------------------------------------------------


# -------------------------- FIT PCF -----------------------
# fit the theoretical PCF parameters to the simulation PCF

fitPCF <- function(){

	# calculate the empirical PCF
	df_PCF <- calcPCF()
	# round distances to integer values
	df_PCF <- df_PCF %>% dplyr::mutate(distance_red=round(distance))
	# calculate mean PCF for reduced distances
	df_PCF_red <- df_PCF %>% dplyr::group_by(distance_red) %>% dplyr::summarise(PCF_mean=mean(PCF))

	# square difference between theory and simulation
	lssqPCF <- function(par,df=df_PCF_red){
		df <- df %>% mutate(PCF_par=PCFparam(distance_red,par)) %>% mutate(PCF_diff=(PCF_mean-PCF_par)^2)
		lssq <- df %>% pull(PCF_diff) %>% sum()
		return(lssq)
	}

	# starting parameter values
	rho <- 100
	lambda <- 10

	# optimise the parameter values -> least square estimate
	optim_par <- optim(par=c(rho,lambda),fn=lssqPCF)

	df_out <- df_PCF_red %>% dplyr::mutate(rho=optim_par$par[1],lambd=optim_par$par[2])

	return(df_out)

}

# ----------------------------------------------------------


# ------------------------ MAKE PLOTS-----------------------

plotSimulation <- function(){

	createDirectory()

	# get the replicate statistics
	df_plot <- readReplicatesDiversity()
	# get grid dimensions
	Nlat <- df_plot$lat %>% unique() %>% length()
	Nlon <- df_plot$lon %>% unique() %>% length()

	# get simulation parameters
	PARs <- readPARAMs()
	hab_switch <- PARs[[2]]
	temp_switch <- PARs[[3]]
	spec_switch <- PARs[[4]]

	# plot the simulation mean
	p <- ggplot() + geom_tile(data=df_plot,aes(x=lon,y=lat,fill=simulation_mean),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='simulation mean',limits=c(0,1)) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10),panel.grid=element_blank()) +
		scale_x_continuous(labels=1:Nlon,breaks=1:Nlon) + scale_y_continuous(labels=1:Nlat,breaks=1:Nlat) 

	p %>% ggsave('Output/plots/simulation_mean.png',.,device='png',width=16,height=15,units='cm')

	# plot the simulation standard deviation
	p <- ggplot() + geom_tile(data=df_plot,aes(x=lon,y=lat,fill=simulation_sd),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='simulation sd') +
		theme(axis.title=element_blank(), legend.title=element_text(size=10),panel.grid=element_blank()) +
		scale_x_continuous(labels=1:Nlon,breaks=1:Nlon) + scale_y_continuous(labels=1:Nlat,breaks=1:Nlat) 

	p %>% ggsave('Output/plots/simulation_sd.png',.,device='png',width=16,height=15,units='cm')

	if(hab_switch){
		# plot the simulation residuals map for habitat area
		p <- ggplot() + geom_tile(data=df_plot,aes(x=lon,y=lat,fill=residualGrad),color='black') + theme_bw() +
			scale_fill_gradient2(low=rgb(1,0,0,0.65),high=rgb(0,1,0,0.65),mid=rgb(1,1,1,0.65),
									name='residuals',limits=c(-1,1),na.value=rgb(0.5,0.5,0.5,0.65)) +
			theme(axis.title=element_blank(), legend.title=element_text(size=10),panel.grid=element_blank()) +
			scale_x_continuous(labels=1:Nlon,breaks=1:Nlon) + scale_y_continuous(labels=1:Nlat,breaks=1:Nlat) 

		p %>% ggsave('Output/plots/residual_grad_map.png',.,device='png',width=16,height=15,units='cm')
	}

	if(temp_switch){
		# plot the simulation residuals map for temperature
		p <- ggplot() + geom_tile(data=df_plot,aes(x=lon,y=lat,fill=residualTemp),color='black') + theme_bw() +
			scale_fill_gradient2(low=rgb(1,0,0,0.65),high=rgb(0,1,0,0.65),mid=rgb(1,1,1,0.65),
									name='residuals',limits=c(-1,1),na.value=rgb(0.5,0.5,0.5,0.65)) +
			theme(axis.title=element_blank(), legend.title=element_text(size=10),panel.grid=element_blank()) +
			scale_x_continuous(labels=1:Nlon,breaks=1:Nlon) + scale_y_continuous(labels=1:Nlat,breaks=1:Nlat) 

		p %>% ggsave('Output/plots/residual_temp_map.png',.,device='png',width=16,height=15,units='cm')
	}

	if(hab_switch){
		# plot the simulation residual distribution for habitat area
		p <- df_plot %>% ggplot() + geom_histogram(aes(residualGrad)) + theme_bw() + xlab('habitat area residuals')
		p %>% ggsave('Output/plots/residual_grad_distribution.png',.,device='png',width=10,height=10,units='cm')
	}

	if(temp_switch){
		# plot the simulation residual distribution for temperature
		p <- df_plot %>% ggplot() + geom_histogram(aes(residualGrad)) + theme_bw() + xlab('temperature residuals')
		p %>% ggsave('Output/plots/residual_temp_distribution.png',.,device='png',width=10,height=10,units='cm')
	}

	if(hab_switch){
		# calculate the regression coefficients between simulation and empirical habitat area
		model <- lm(data=df_plot,simulation_mean~gradNorm)
		coeff <- summary(model)$coefficients %>% as.data.frame() %>% dplyr::select(Estimate)
		# calculate the Pearson correlation coefficient
		correl <- cor.test(df_plot$gradNorm,df_plot$simulation_mean,method='pearson')$estimate %>% as.numeric()
		# plot the scatter
		p <- df_plot %>% ggplot() + geom_point(aes(x=gradNorm,y=simulation_mean),pch=19) + 
						geom_abline(intercept=coeff[1,1], slope=coeff[2,1], color='red', linewidth=1) +
						geom_abline(intercept=0,slope=1,linetype=2,linewidth=0.75) +
						xlab('habitat gradient') + ylab('simulation') +
						theme_bw() + xlim(0,1) + ylim(0,1) + labs(title=paste0('r = ',round(correl,4)))
		p %>% ggsave('Output/plots/scatter_grad.png',.,device='png',width=10,height=10,units='cm')
	}

	if(temp_switch){
		# calculate the regression coefficients between simulation and empirical temperature
		model <- lm(data=df_plot,simulation_mean~boltzNorm)
		coeff <- summary(model)$coefficients %>% as.data.frame() %>% dplyr::select(Estimate)
		# calculate the Pearson correlation coefficient
		correl <- cor.test(df_plot$boltzNorm,df_plot$simulation_mean,method='pearson')$estimate %>% as.numeric()
		# plot the scatter
		p <- df_plot %>% ggplot() + geom_point(aes(x=boltzNorm,y=simulation_mean),pch=19) + 
						geom_abline(intercept=coeff[1,1], slope=coeff[2,1], color='red', linewidth=1) +
						geom_abline(intercept=0,slope=1,linetype=2,linewidth=0.75) +
						xlab('Boltzmann factor') + ylab('simulation') +
						theme_bw() + xlim(0,1) + ylim(0,1) 
		p %>% ggsave('Output/plots/scatter_temp.png',.,device='png',width=10,height=10,units='cm')

	}

	if(spec_switch){
		# get the fitted PCF
		df_PCF_red <- fitPCF()
		Rho <- df_PCF_red$rho %>% mean()
		Lambd <- df_PCF_red$lambd %>% mean()
		# dummy distances
		x <- seq(min(df_PCF_red$distance_red),max(df_PCF_red$distance_red),0.01)
		# calculate theoretical prediction
		df_PCF_fit <- data.frame(x=x,y=PCFparam(x,c(Rho,Lambd)))
		# plot the reduced PCF
		p <- df_PCF_red %>% ggplot() + geom_point(aes(x=distance_red,y=PCF_mean)) + 
									geom_line(data=df_PCF_fit,aes(x=x,y=y)) +
									theme_bw() + xlab('distance') + ylab('PCF') + labs(title=paste0('ρ = ',round(Rho,4),'   ','λ = ',round(Lambd,4))) +
									ylim(0,max(c(df_PCF_red$PCF_mean,df_PCF_fit$y)))

		p %>% ggsave('Output/plots/PCF_reduced.png',.,device='png',width=15,height=10,units='cm')
	}

}

# ----------------------------------------------------------


# ---------------- MAKE THE PLOTS --------------------------
t1 <- Sys.time()

plotSimulation()

t2 <- Sys.time()

td <- seconds_to_period(difftime(t2,t1,units='secs'))
print(paste0('wall time: ',sprintf('%02d:%02d:%02d',td@hour, minute(td), round(second(td))) ))

# ----------------------------------------------------------


# --------------- CREATE THE SIMULATION GRID ---------------
# Nlon = number of longitudinal cells
# Nlat = number of latitudinal cells
# gradIn = max(habitat) - min(habitat)
# gradTemp = max(temperature) - min(temperature)

createGrid <- function(Nlon,Nlat,gradIn,gradTemp){

	# create the coordinate lists
	lon_list <- rep(1:Nlon,Nlat)
	lat_list <- c()
	for(i in 1:Nlat){
		lat_list <- c(lat_list,rep(i,Nlon))
	}

	# difference in habitat area between neighbouring rings:
	# . . . . .
	# . x x x .
	# . x o x .
	# . x x x .
	# . . . . .
	delta_grad <- (gradIn-1) / ((Nlon-1)/2 - 1)

	# create grid container
	df_out <- data.frame(lat=lat_list,lon=lon_list,grad=NA,temp=NA)

	# for each ring around centre, starting from centre (i = 0)
	for(i in 0:(gradIn-2)){
		# value of the habitat area
		grad_i <- gradIn - i*delta_grad
		# create the concentric rings
		df_out$grad[abs(df_out$lon - (Nlon+1)/2) == i & abs(df_out$lat - (Nlat+1)/2) %in% 0:i] <- grad_i
		df_out$grad[abs(df_out$lon - (Nlon+1)/2) %in% 0:i & abs(df_out$lat - (Nlat+1)/2) == i ] <- grad_i
	}

	# create the latitudinal temperature gradient:
	# . . .
	# x x x
	# o o o
	for(i in 1:Nlat){
		df_out$temp[df_out$lat==i] <- gradTemp/(1-Nlat)*(i-Nlat)
	}

	# remaining habitat area = 1
	df_out$grad[is.na(df_out$grad)] <- 1
	# transform temperature to Kelvin
	df_out$temp <- df_out$temp + 273.15
	
	return(df_out)

}

# ----------------------------------------------------------


# ------------- CREATE THE RANDOM SIMULATION GRID ----------
# Nlon = number of longitudinal cells
# Nlat = number of latitudinal cells
# gradIn = max(habitat) - min(habitat)
# gradTemp = max(temperature) - min(temperature)

createGridRandom <- function(Nlon,Nlat,gradIn,gradTemp){

	# create coordinate list
	lon_list <- rep(1:Nlon,Nlat)
	lat_list <- c()
	for(i in 1:Nlat){
		lat_list <- c(lat_list,rep(i,Nlon))
	}

	# random habitat area and temperature values
	df_out <- data.frame(lat=lat_list,lon=lon_list,
							grad=runif(length(lat_list),min=1,max=gradIn),
							temp=runif(length(lat_list),min=0,max=gradTemp))	

	return(df_out)
}

# ----------------------------------------------------------

# create examples of simulation grids:
# uncomment to write out the grid

#write.csv(createGrid(11,11,5,30),'GridFiles/grid_test_scaling.csv')
#write.csv(createGridRandom(11,11,5,30),'GridFiles/grid_test_scaling_rand.csv')

# ---------------- PLOT SIMULATION GRID --------------------
# plot all simulation grids in GridFiles

plotGrid <- function(){

	# list of all simulation grids in GridFiles
	list_files <- list.files(path='GridFiles',pattern='*.csv',full.names=TRUE)

	for(file in list_files){

		# get the grid
		df <- read.csv(file,header=TRUE,row.names=1)

		# get grid specifics
		Nlon <- df %>% dplyr::pull(lon) %>% unique() %>% length()
		Nlat <- df %>% dplyr::pull(lat) %>% unique() %>% length()
		gradIn <- df %>% dplyr::pull(grad) %>% max()

		grid_name <- str_split(file,pattern='/')[[1]][2]
		grid_name <- str_split(grid_name,pattern='[.]')[[1]][1]

		p_hab <- ggplot() + geom_tile(data=df,aes(x=lon,y=lat,fill=grad),color='black') + theme_bw() +
					scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='habitat area',limits=c(0,gradIn)) +
					scale_x_continuous(labels=1:Nlon,breaks=1:Nlon) + scale_y_continuous(labels=1:Nlat,breaks=1:Nlat) +
					theme(axis.title=element_blank(), panel.grid=element_blank())

		p_temp <- ggplot() + geom_tile(data=df,aes(x=lon,y=lat,fill=temp),color='black') + theme_bw() +
					scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='temperature (K)') +
					scale_x_continuous(labels=1:Nlon,breaks=1:Nlon) + scale_y_continuous(labels=1:Nlat,breaks=1:Nlat) +
					theme(axis.title=element_blank(),panel.grid=element_blank())

		# Create a text grob
		tgrob <- text_grob(grid_name,size = 14)
		# Draw the main title = grid name
		p_text <- as_ggplot(tgrob) + theme(plot.margin = margin(0,5,0,0, "cm"))

		q <- ggarrange(p_text,NULL,p_hab,p_temp,nrow=2,ncol=2,heights=c(1,6))

		q %>% ggsave(paste0('GridFiles/',grid_name,'.png'),.,device='png',width=24,height=10,units='cm',bg='white')

	}

}

# uncomment to plot the simulation grids
#plotGrid()

# ----------------------------------------------------------