# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Simulation of Mechanistic Model
#####################################
# MAIN EXECUTION of the simulation
#####################################


# ------------ IMPORT FROM OTHER FILES ---------------------

from Grid import*
from parameters import*

# ----------------------------------------------------------


# ----------------- IMPORT MODULES -------------------------

import time
import os
import shutil

# ----------------------------------------------------------

# create the simulation grid
g = Grid()
# fill the grid 
# initial number of species Nspec0 defined in parameters.py
g.fillGrid(Nspec0)
# print the initial grid before the first turnover
g.printGrid(0)

# create a directory for output files
if(os.path.exists("Output")):
	shutil.rmtree("Output")
os.mkdir("Output")

# time the complete run
t1 = time.time()
# for each turnover
# number of turnovers in simulation Nturn defined in paramaters.py
for i in range(0,Nturn):
	print('turnover ' + str(i) + ' of ' + str(Nturn),end='\r')
	# turn over the grid
	g.turnover()

	# every Nout turnovers:
	# Nout defined in parameters.py
	if(i%Nout==0):
		# update the species list
		g.updateSpecies()
		# print the simulation grid
		g.printGrid(i/Nout)

t2 = time.time()
print('wall time: ' + str(round((t2-t1)/60,4)) + ' minutes')
print_parameters(t2-t1)

# ----------------------------------------------------------