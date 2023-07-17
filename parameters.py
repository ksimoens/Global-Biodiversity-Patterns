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

# number of cells
# collected directly from the Worm and Tittensor (2018) file
Nlon = len(np.unique(pd.read_csv('land_all_data_inc_working.csv').iloc[:,[5]])) # longitudinal
Nlat = len(np.unique(pd.read_csv('land_all_data_inc_working.csv').iloc[:,[6]])) # latitudinal
# number of populations in a grid cell
Nloc = 16
# migration rate
Pdisp = 0.1
# speciation rate
Pspec = 0.01
# initial number of species
Nspec0 = 1
# number of turnovers in the simulation
Nturn = 5000000
# number of turnovers after which to print the grid
Nout = 5000

# ----------------------------------------------------------

# print the simulation specifics
def print_parameters(runtime):

	with open('Output/PARAM_file.txt', 'w') as f:
		f.write('####################\n')
		f.write('--BASIC SIMULATION--\n')
		f.write('####################\n')
		now = datetime.datetime.now()
		f.write(now.strftime("%d/%m/%Y %H:%M:%S") + '\n')
		f.write('####################\n')
		f.write('\n')
		f.write('Main grid characteristics:\n')
		f.write('\t- global grid\n')
		f.write('\t- periodic longitudinal boundary conditions\n')
		f.write('\t- hard latitudinal boundary conditions\n')
		f.write('\t- purely neutral dynamics\n')
		f.write('\n')
		f.write('####################\n')
		f.write('\n')
		f.write('Simulation parameters:\n')
		f.write('\n')
		f.write('Number of longitudinal cells: ' + str(Nlon) + '\n')
		f.write('Number of latitudinal cells: ' + str(Nlat) + '\n')
		f.write('Number of local populations: ' + str(Nloc) + '\n')
		f.write('Migration rate: ' + str(Pdisp) + '\n')
		f.write('Speciation rate: ' + str(Pspec) + '\n')
		f.write('Initial number of species: ' + str(Nspec0) + '\n')
		f.write('Number of simulated turnovers: ' + str(Nturn) + '\n')
		f.write('Number of turnovers after which to output: ' + str(Nout) + '\n')
		f.write('\n')
		f.write('####################\n')
		f.write('\n')
		f.write('Total runtime: ' + str(datetime.timedelta(seconds=runtime)) + '\n')
		f.write('\n')
		f.write('####################\n')
		print(runtime)

		f.close()