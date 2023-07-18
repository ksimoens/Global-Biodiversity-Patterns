# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Simulation of Mechanistic Model
#####################################
# PARAMETERS of the simulation
#####################################


# ----------------- IMPORT MODULES -------------------------

import pandas as pd
import numpy as np
import datetime

# ----------------------------------------------------------

# define the name of the simulation grid file
# this file is located in the GridFiles directory
GridFile = 'CPR_grid_env_div_sim.csv'

# number of cells
# collected directly from the simulation grid file
Nlon = len(np.unique(pd.read_csv('GridFiles/' + GridFile)['x'])) # longitudinal
Nlat = len(np.unique(pd.read_csv('GridFiles/' + GridFile)['y'])) # latitudinal
# number of populations in a grid cell
Nloc = 16
# migration rate
Pdisp = pow(10,-1.5)
# speciation rate
Pspec = 0.1
# width of the species temperature range
NicheWidth = 10
# slope of the linear relationship between the habitat area proxy and the local number of populations
# if the habitat area proxy x 2 => local number of populations x HabSlope
HabSlope = 2
# number of computer cores to use in parallellisation
NCPU = 4
# number of replicate runs to disperse over the computer cores
Nrep = 20

# switches to implement the new physics
# True = active
# False = inactive
TempTurnover = True # temperature-dependent death rates
TempSpeciation = True # temperature-dependent speciation rates
AreaHabitat = False # scaling of local number of populations with habitat area
TempNiches = False # temperature ranges and migration limitation

# print species information in the output
# True = print full diversity matrix
# False = print only a local diversity column
# not printing the diversity matrix is a lot faster
SpeciesInformation = False


# ----------------------------------------------------------


# print the simulation specifics
def print_parameters(runtime):

	with open('Output/PARAM_file.txt', 'w') as f:
		f.write('####################\n')
		f.write('--DATA SIMULATION--\n')
		f.write('####################\n')
		now = datetime.datetime.now()
		f.write(now.strftime("%d/%m/%Y %H:%M:%S") + '\n')
		f.write('####################\n')
		f.write('\n')
		f.write('Main grid characteristics:\n')
		f.write('\t- grid: ' + GridFile + '\n')
		f.write('\t- hard longitudinal boundary conditions\n')
		f.write('\t- hard latitudinal boundary conditions\n')
		f.write('\n')
		f.write('####################\n')
		f.write('\n')
		f.write('Numerical specifications:\n')
		f.write('\n')
		f.write('\t- Coalescence algorithm\n')
		f.write('\t- Parallel computing\n')
		f.write('\t- Number of computing cores: ' + str(NCPU) + '\n')
		f.write('\t- Number of replicate runs: ' + str(Nrep) + '\n')
		f.write('\n')
		f.write('####################\n')
		f.write('\n')
		f.write('Simulation parameters:\n')
		f.write('\n')
		f.write('Number of longitudinal cells: ' + str(Nlon) + '\n')
		f.write('Number of latitudinal cells: ' + str(Nlat) + '\n')
		f.write('Base number of local populations: ' + str(Nloc) + '\n')
		f.write('Migration rate: ' + str(Pdisp) + '\n')
		f.write('Base speciation rate: ' + str(Pspec) + '\n')
		if(TempTurnover):
			f.write('Temperature-dependent death rates: active\n')
		else:
			f.write('Temperature-dependent death rates: inactive\n')
		if(TempSpeciation):
			f.write('Temperature-dependent speciation rates: active\n')
		else:
			f.write('Temperature-dependent speciation rates: inactive\n')
		if(AreaHabitat):
			f.write('Habitat area scaling: active\n')
			f.write('Slope of habitat scaling relationship: ' + str(HabSlope) + '\n')
		else:
			f.write('habitat area scaling: inactive\n')
		if(TempNiches):
			f.write('temperature ranges of species: active\n')
		else:
			f.write('temperature ranges of species: inactive\n')
		if(SpeciesInformation):
			f.write('print complete diversity matrix to output: active\n')
		else:
			f.write('print complete diversity matrix to output: inactive\n')
		f.write('\n')
		f.write('####################\n')
		f.write('\n')
		f.write('Total runtime: ' + str(datetime.timedelta(seconds=runtime)) + '\n')
		f.write('\n')
		f.write('####################\n')
		print(runtime)

		f.close()