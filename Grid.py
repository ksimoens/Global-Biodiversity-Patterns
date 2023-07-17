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
		
		# import the grid cell coordinates
		# 'land_all_data_inc_working.csv' from Worm and Tittensor (2018)
		df_grid = pd.read_csv('land_all_data_inc_working.csv').iloc[:,[0,5,6]]
		# fill the grid with Local objects
		# Nlon / Nlat defined in parameters.py
		k = 0
		for i in range(0,Nlon):
			for j in range(0,Nlat):
				self.global_grid.append( Local(df_grid.iloc[k,1],df_grid.iloc[k,2],k+1) )
				k += 1

		# total number of species in the grid
		self.Nspec = 0
		# list of species in the grid: list of integers
		self.species = np.zeros(self.Nspec)
		# maximum integer value in the species list
		self.MaxSpec = 0

	# fill the grid befor the simulation
	# requires an initial number of species: Nspec
	def fillGrid(self,Nspec):
		
		# fill each Local object with populations
		# fillLocal defined in Local.py
		for loc in self.global_grid:
			loc.fillLocal(Nspec)
		# update species list
		self.Nspec = Nspec
		self.species = np.arange(0,Nspec)
		self.MaxSpec = np.max(self.species)

	# after a few turnovers, update the species list
	def updateSpecies(self):

		# list species identity of all local populations
		pop_list = np.array([])
		for cell in self.global_grid:
			pop_list = np.concatenate((pop_list,cell.populations))

		# list the identity of each population
		# add present species to the new species list
		# update species information
		self.species = np.array([])
		for i in range(0,self.MaxSpec+1):
			n = np.count_nonzero(pop_list == i)
			if(n != 0):
				self.species = np.append(self.species,i)

		self.Nspec = len(self.species)

	# print grid information to a csv file
	# requires an identifier 'step' to name the output file
	# best to use updateSpecies() before print
	def printGrid(self,step):

		# initialise lat and lon containers
		lon_list = np.zeros(len(self.global_grid))
		lat_list = np.zeros(len(self.global_grid))
		# initialise species container
		spec_list = np.zeros( (len(self.global_grid),self.Nspec) )

		# for each grid cell, collect diversity matrix
		# can take a while for many species
		for i in range(0,len(lat_list)):
			lon_list[i] = self.global_grid[i].lon
			lat_list[i] = self.global_grid[i].lat

			for j in range(0,self.Nspec):
				spec_list[i][j] = np.count_nonzero(self.global_grid[i].populations == self.species[j])

		# create pandas dataframe and write out to csv file
		darray = np.concatenate((lon_list, lat_list)).reshape((-1, 2), order='F')
		darray = np.concatenate( (darray,spec_list) ,axis=1)
		names_out = ['lon','lat'] + ['spec_' + str(int(ID)) for ID in self.species] 		
		df_out = pd.DataFrame(data=darray, columns=names_out)
		
		df_out.to_csv('Output/grid_' + str(int(step)).zfill(4) + '.csv')

	# get Moore neighbourhood of cell with index
	# assumes periodic longitudinal boundary conditions 
	# assumes hard latitudinal boundary conditions
	# returns list of integers representing the species identity of local populations
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

		if(index >= Nlon):
			# 3
			# - x -
			# - o -
			# - - -
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon].populations))
			
			# 4
			# | x - . - - |
			# | - - . - o |
			# | - - . - - |
			if((index+1) % Nlon == 0):
				neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon-Nlon+1].populations))
			# - - x
			# - o -
			# - - -
			else:
				neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon+1].populations))

			# 5
			# | - - . - x |
			# | o - . - - |
			# | - - . - - |
			if((index) % 45 == 0):
				neighbours = np.concatenate((neighbours , self.global_grid[index-1].populations))
			# x - -
			# - o -
			# - - -
			else:
				neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon-1].populations))
	#	else: nothing
			# _ _ _
			# - o -
			# - - -

		if(index < len(self.global_grid)-Nlon):
			# 6
			# - - -
			# - o -
			# - x -
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon].populations))

			# 7
			# | - - . - - |
			# | - - . - o |
			# | x - . - - |
			if((index+1) % Nlon == 0):
				neighbours = np.concatenate((neighbours , self.global_grid[index+1].populations))
			# - - -
			# - o -
			# - - x
			else:
				neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon+1].populations))

			# 8
			# | - - . - - |
			# | o - . - - |
			# | - - . - x |
			if((index) % 45 == 0):
				neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon+Nlon-1].populations))
			# - - - 
			# - o -
			# x - -
			else:
				neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon-1].populations))
			
	#	else: nothing
			# - - -
			# - o -
			# _ _ _

		return(neighbours)

	# select a random cell and return index
	def selectCell(self):

		i = random.randint(0,len(self.global_grid)-1)

		return(i)

	# do a single forward turnover
	def turnover(self):

		# select a random cell to be disturbed
		i = self.selectCell()

		# select local population to die
		old = random.randint(0,len(self.global_grid[i].populations)-1)

		# random number to check for migration
		r = random.uniform(0,1)
		# initialise pool of possible dispersers
		disp_pool = np.array([])

		# check for migration
		# Pdisp defined in parameters.py
		if(r < Pdisp):
			# migration: pool of dispersers = Moore neighbourhood
			disp_pool = self.getNeighbours(i)
		else:
			# no migration: pool of dispersors = local populations
			disp_pool = np.delete(self.global_grid[i].populations,old)

		# random number to check for speciation
		s = random.uniform(0,1)

		# check for speciation
		# Pspec defined in parameters.py
		if(s < Pspec):
			# update the species list
			self.global_grid[i].populations[old] = self.MaxSpec + 1
			self.MaxSpec += 1
		else:
			# select a disperser to take the free spot
			new = random.randint(0,len(disp_pool)-1)
			self.global_grid[i].populations[old] = disp_pool[new]

# ----------------------------------------------------------