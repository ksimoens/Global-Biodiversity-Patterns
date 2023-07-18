# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Simulation of Mechanistic Model
#####################################
# SIMULATION GRID and related methods
#####################################

# ------------ IMPORT FROM OTHER FILES ---------------------

from parameters import*
from Local import*

# ----------------------------------------------------------


# ----------------- IMPORT MODULES -------------------------

import numpy as np
import pandas as pd
import random
import os
import shutil

# ----------------------------------------------------------


# ---------------------- GRID CLASS ------------------------
class Grid():

	# define constructor	
	def __init__(self):

		# global grid = list of Local objects
		# Local defined in Local.py
		self.global_grid = []
		
		# import the constructed simulation grid
		# constructed in https://github.com/ksimoens/Thesis-Data-Analysis.git
		# GridFile defined in parameters.py
		df_grid = pd.read_csv('GridFiles/' + GridFile)

		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		# extract/calculate the habitat area proxy
		# for now this has to be done explicitly here
		# a better solution is required
		# AreaHabitat defined in parameters.py
		if(AreaHabitat):
			# example: combine forest biomass and altitudinal gradient and take 10 base logarithm
			foreGrad = np.asarray(df_grid['fore']*df_grid['grad'])
			foreGrad = foreGrad / np.nanmin(foreGrad) 
			hab = np.log10(foreGrad) + 1.
		else:
			# just use an empty array if habitat area is not used
			hab = np.empty(len(df_grid))
		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		# add all relevant columns in order to define the simulation grid
		# each Local object remembers all this information during the simulation
		k = 0
		for i in range(0,Nlon):
			for j in range(0,Nlat):
				self.global_grid.append( Local(df_grid.loc[k, 'x'],df_grid.loc[k, 'y'],k,df_grid.loc[k, 'active'],df_grid.loc[k, 'upper'],df_grid.loc[k, 'lower'],df_grid.loc[k, 'left'],df_grid.loc[k, 'right'],df_grid.loc[k, 'upper_left'],df_grid.loc[k, 'upper_right'],df_grid.loc[k, 'lower_left'],df_grid.loc[k, 'lower_right'],df_grid.loc[k, 'temp'],hab[k]) )
				k += 1

		# lowest temperature in the grid
		self.temp_min = df_grid.loc[:,'temp'].abs().min()

	# fill the grid befor the simulation
	def fillGrid(self):
		
		# fill each Local object with populations
		# fillLocal defined in Local.py
		for loc in self.global_grid:
			if(loc.active == 1):
				loc.fillLocal()

	# get Moore neighbourhood of cell with index
	# assumes hard longitudinal boundary conditions 
	# assumes hard latitudinal boundary conditions
	# returns list of Population objects 
	# Nlon defined in parameters.py
	def getNeighbours(self,index):

		neighbours = np.array([])

		# get the geometrical booleans 
		# convey whether a neighbour is active
		up = self.global_grid[index].upper
		low = self.global_grid[index].lower
		lef = self.global_grid[index].left
		rig = self.global_grid[index].right
		up_lef = self.global_grid[index].upper_left
		up_rig = self.global_grid[index].upper_right
		low_lef = self.global_grid[index].lower_left
		low_rig = self.global_grid[index].lower_right

		# 1
		# - - -
		# - o x
		# - - -
		if(rig == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index+1].populations))

		# 2
		# - - -
		# x o -
		# - - -
		if(lef == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index-1].populations))
		
		# 3
		# - x -
		# - o -
		# - - -
		if(up == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon].populations))
		
		# 4	
		# - - x
		# - o -
		# - - -
		if(up_rig == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon+1].populations))

		# 5	
		# x - -
		# - o -
		# - - -
		if(up_lef == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon-1].populations))
	
		# 6
		# - - -
		# - o -
		# - x -
		if(low == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon].populations))

		# 7
		# - - -
		# - o -
		# - - x
		if(low_rig == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon+1].populations))

		# 8
		# - - - 
		# - o -
		# x - -
		if(low_lef == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon-1].populations))
			

		return(neighbours)