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
		# constructed in plot.R
		# GridFile defined in parameters.py		
		df_grid = pd.read_csv('GridFiles/' + GridFile,index_col=0)

		# add all relevant columns in order to define the simulation grid
		# each Local object remembers all this information during the simulation
		k = 0
		for i in range(0,Nlon):
			for j in range(0,Nlat):
				self.global_grid.append( Local(df_grid.iloc[k,1],df_grid.iloc[k,0],k,df_grid.iloc[k,2],df_grid.iloc[k,3]) ) 
				k += 1

		# lowest temperature in the grid
		self.Tmin = df_grid['temp'].min()

	# fill the grid befor the simulation
	def fillGrid(self):
		
		# fill each Local object with populations
		# fillLocal defined in Local.py
		for loc in self.global_grid:
			loc.fillLocal()

	
	# get Moore neighbourhood of cell with index
	# assumes periodical longitudinal boundary conditions 
	# assumes periodical latitudinal boundary conditions
	# returns list of Population objects 
	# Nlon/Nlat defined in parameters.py
	def getNeighbours(self,index):

		neighbours = np.array([])

		# 1
		# | - - . - - |
		# | x - . - o |
		# | - - . - - |
		if((index+1) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon+1].populations))
		# - - -
		# - o x
		# - - -
		else:
			neighbours = np.concatenate((neighbours , self.global_grid[index+1].populations))

		# 2
		# | - - . - - |
		# | o - . - x |
		# | - - . - - |
		if((index) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon-1].populations))
		# - - -
		# x o -
		# - - -
		else:
			neighbours = np.concatenate((neighbours , self.global_grid[index-1].populations))

		# 3
		# _____
		# - o -
		# - - -
		# . . .
		# - - -
		# - x -
		# _____
		if(index < Nlon):
			neighbours = np.concatenate((neighbours , self.global_grid[Nlon*Nlat-Nlon+index].populations))
		# - x -
		# - o -
		# - - -
		else:
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon].populations))

		# 4
		#  __________
		# | - - . - o |
		# | - - . - - |
		# | . . . . . |
		# | - - . - - |
		# | x - . - - |
		#  ___________
		if(index < Nlon and (index+1) % Nlon == 0):	
			neighbours = np.concatenate((neighbours , self.global_grid[Nlon*Nlat-Nlon].populations))
		# | x - . - - |
		# | - - . - o |
		# | - - . - - |
		elif(index >= Nlon and (index+1) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon-Nlon+1].populations))
		# _____
		# - o -
		# - - - 
		# . . .
		# - - - 
		# - - x
		# _____
		elif(index < Nlon and (index+1) % Nlon > 0):
			neighbours = np.concatenate((neighbours , self.global_grid[Nlon*Nlat-Nlon+index+1].populations))
		# - - x
		# - o -
		# - - -
		else:
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon+1].populations))

		# 5
		#  ___________
		# | o - . - - |
		# | - - . - - |
		# | . . . . . |
		# | - - . - - |
		# | - - . - x |
		#  ___________
		if(index < Nlon and (index) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[Nlon*Nlat-1].populations))

		# | - - . - x |
		# | o - . - - |
		# | - - . - - |
		elif(index >= Nlon and (index) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index-1].populations))
		# _____
		# - o -
		# - - -
		# . . .
		# - - -
		# x - -
		# _____
		elif(index < Nlon and (index) % Nlon > 0):
			neighbours = np.concatenate((neighbours , self.global_grid[Nlon*Nlat-1-Nlon+index].populations))
		# x - -
		# - o -
		# - - -
		else:
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon-1].populations))
	
		# 6
		# _____
		# - x -
		# - - -
		# . . .
		# - - -
		# - o -
		# _____
		if(index >= len(self.global_grid)-Nlon):
			neighbours = np.concatenate((neighbours , self.global_grid[index-(Nlon*Nlat-Nlon)].populations))
		# - - -
		# - o -
		# - x -
		else:
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon].populations))

		# 7
		#  ___________
		# | x - . - - |
		# | - - . - - |
		# | . . . . . |
		# | - - . - - |
		# | - - . - o |
		#  ___________
		if(index >= len(self.global_grid)-Nlon and (index+1) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index - (Nlon*Nlat-1)].populations))
		# | - - . - - |
		# | - - . - o |
		# | x - . - - |
		elif(index < len(self.global_grid)-Nlon and (index+1) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index+1].populations))
		# _____
		# - - x
		# - - -
		# . . .
		# - - -
		# - o -
		# _____
		elif(index >= len(self.global_grid)-Nlon and (index+1) % Nlon > 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index-(Nlon*Nlat-Nlon)+1].populations))
		# - - -
		# - o -
		# - - x
		else:
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon+1].populations))

		# 8
		#  ___________
		# | - - . - x |
		# | - - . - - |
		# | . . . . . |
		# | - - . - - |
		# | o - . - - |
		#  ___________
		if(index >= len(self.global_grid)-Nlon and (index) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[Nlon - 1 - (index - (Nlon*Nlat-Nlon) )].populations))
		# | - - . - - |
		# | o - . - - |
		# | - - . - x |
		elif(index < len(self.global_grid)-Nlon and (index) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon+Nlon-1].populations))
		# _____
		# x - -
		# - - -
		# . . .
		# - - -
		# - o -
		# _____
		elif(index >= len(self.global_grid)-Nlon and (index) % Nlon > 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index-(Nlon*Nlat-Nlon)-1].populations))
		# - - - 
		# - o -
		# x - -
		else:
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon-1].populations))

		return(neighbours)


	# same as in getNeighbours, but now return a list of neighbour indices
	# used in an alternative migration simulation
	def getNeighboursIndex(self,index):

		neighbours = np.array([])

		# 1
		# | - - . - - |
		# | x - . - o |
		# | - - . - - |
		if((index+1) % Nlon == 0):
			neighbours = np.append(neighbours,index-Nlon+1)
		# - - -
		# - o x
		# - - -
		else:
			neighbours = np.append(neighbours,index+1)

		# 2
		# | - - . - - |
		# | o - . - x |
		# | - - . - - |
		if((index) % Nlon == 0):
			neighbours = np.append(neighbours,index+Nlon-1)
		# - - -
		# x o -
		# - - -
		else:
			neighbours = np.append(neighbours,index-1)

		# 3
		# _____
		# - o -
		# - - -
		# . . .
		# - - -
		# - x -
		# _____
		if(index < Nlon):
			neighbours = np.append(neighbours, Nlon*Nlat-Nlon+index)
		# - x -
		# - o -
		# - - -
		else:
			neighbours = np.append(neighbours,index-Nlon)

		# 4
		#  __________
		# | - - . - o |
		# | - - . - - |
		# | . . . . . |
		# | - - . - - |
		# | x - . - - |
		#  ___________
		if(index < Nlon and (index+1) % Nlon == 0):	
			neighbours = np.append(neighbours,Nlon*Nlat-Nlon)
		# | x - . - - |
		# | - - . - o |
		# | - - . - - |
		elif(index >= Nlon and (index+1) % Nlon == 0):
			neighbours = np.append(neighbours,index-Nlon-Nlon+1)
		# _____
		# - o -
		# - - - 
		# . . .
		# - - - 
		# - - x
		# _____
		elif(index < Nlon and (index+1) % Nlon > 0):
			neighbours = np.append(neighbours,Nlon*Nlat-Nlon+index+1)
		# - - x
		# - o -
		# - - -
		else:
			neighbours = np.append(neighbours,index-Nlon+1)

		# 5
		#  ___________
		# | o - . - - |
		# | - - . - - |
		# | . . . . . |
		# | - - . - - |
		# | - - . - x |
		#  ___________
		if(index < Nlon and (index) % Nlon == 0):
			neighbours = np.append(neighbours,Nlon*Nlat-1)
		# | - - . - x |
		# | o - . - - |
		# | - - . - - |
		elif(index >= Nlon and (index) % Nlon == 0):
			neighbours = np.append(neighbours,index-1)
		# _____
		# - o -
		# - - -
		# . . .
		# - - -
		# x - -
		# _____
		elif(index < Nlon and (index) % Nlon > 0):
			neighbours = np.append(neighbours,Nlon*Nlat-1-Nlon+index)
		# x - -
		# - o -
		# - - -
		else:
			neighbours = np.append(neighbours,index-Nlon-1)
	
		# 6
		# _____
		# - x -
		# - - -
		# . . .
		# - - -
		# - o -
		# _____
		if(index >= len(self.global_grid)-Nlon):
			neighbours = np.append(neighbours,index-(Nlon*Nlat-Nlon))
		# - - -
		# - o -
		# - x -
		else:
			neighbours = np.append(neighbours,index+Nlon)

		# 7
		#  ___________
		# | x - . - - |
		# | - - . - - |
		# | . . . . . |
		# | - - . - - |
		# | - - . - o |
		#  ___________
		if(index >= len(self.global_grid)-Nlon and (index+1) % Nlon == 0):
			neighbours = np.append(neighbours,index - (Nlon*Nlat-1))
		# | - - . - - |
		# | - - . - o |
		# | x - . - - |
		elif(index < len(self.global_grid)-Nlon and (index+1) % Nlon == 0):
			neighbours = np.append(neighbours,index+1)
		# _____
		# - - x
		# - - -
		# . . .
		# - - -
		# - o -
		# _____
		elif(index >= len(self.global_grid)-Nlon and (index+1) % Nlon > 0):
			neighbours = np.append(neighbours,index-(Nlon*Nlat-Nlon)+1)
		# - - -
		# - o -
		# - - x
		else:
			neighbours = np.append(neighbours,index+Nlon+1)

		# 8
		#  ___________
		# | - - . - x |
		# | - - . - - |
		# | . . . . . |
		# | - - . - - |
		# | o - . - - |
		#  ___________
		if(index >= len(self.global_grid)-Nlon and (index) % Nlon == 0):
			neighbours = np.append(neighbours,Nlon - 1 - (index - (Nlon*Nlat-Nlon) ))
		# | - - . - - |
		# | o - . - - |
		# | - - . - x |
		elif(index < len(self.global_grid)-Nlon and (index) % Nlon == 0):
			neighbours = np.append(neighbours,index+Nlon+Nlon-1)
		# _____
		# x - -
		# - - -
		# . . .
		# - - -
		# - o -
		# _____
		elif(index >= len(self.global_grid)-Nlon and (index) % Nlon > 0):
			neighbours = np.append(neighbours,index-(Nlon*Nlat-Nlon)-1)
		# - - - 
		# - o -
		# x - -
		else:
			neighbours = np.append(neighbours,index+Nlon-1)

		neighbours = neighbours.astype(int)

		return(neighbours)

	# return all populations outside the Moore neighbourhood of and index
	# used in an alternative migration simulation 
	def getAllPopulations(self,index):

		# get the indices of the Moore neighbourhood
		indices = np.append(self.getNeighboursIndex(index),index)
		# container for all populations
		all_pops = np.array([])
		# add all populations that are not in the Moore neighbourhood
		for i in range(0,len(self.global_grid)):
			if(i not in indices):
				all_pops = np.concatenate((all_pops,self.global_grid[i].populations))

		return(all_pops)