import numpy as np
import pandas as pd
import random
import time
import os
import shutil

from parameters import*
from Local import*

class Grid():

	def __init__(self):

		self.global_grid = []
		
		df_grid = pd.read_csv('land_all_data_inc_working.csv').iloc[:,[0,5,6]]
		k = 0
		for i in range(0,Nlon):
			for j in range(0,Nlat):
				self.global_grid.append( Local(df_grid.iloc[k,1],df_grid.iloc[k,2],k+1) )
				k += 1

		self.Nspec = 0
		self.species = np.zeros(self.Nspec)
		self.MaxSpec = 0

	def fillGrid(self,Nspec):
		
		for loc in self.global_grid:
			loc.fillLocal(Nspec)
		self.Nspec = Nspec
		self.species = np.arange(0,Nspec)
		self.MaxSpec = np.max(self.species)

	def updateSpecies(self):

		pop_list = np.array([])
		for cell in self.global_grid:
			pop_list = np.concatenate((pop_list,cell.populations))

		self.species = np.array([])
		for i in range(0,self.MaxSpec+1):
			n = np.count_nonzero(pop_list == i)
			if(n != 0):
				self.species = np.append(self.species,i)

		self.Nspec = len(self.species)
		print(self.species)

	def printGrid(self,step):

		lon_list = np.zeros(len(self.global_grid))
		lat_list = np.zeros(len(self.global_grid))

		spec_list = np.zeros( (len(self.global_grid),self.Nspec) )

		for i in range(0,len(lat_list)):
			lon_list[i] = self.global_grid[i].lon
			lat_list[i] = self.global_grid[i].lat

			for j in range(0,self.Nspec):
				spec_list[i][j] = np.count_nonzero(self.global_grid[i].populations == self.species[j])


		darray = np.concatenate((lon_list, lat_list)).reshape((-1, 2), order='F')
		darray = np.concatenate( (darray,spec_list) ,axis=1)
		names_out = ['lon','lat'] + ['spec_' + str(int(ID)) for ID in self.species] 		
		df_out = pd.DataFrame(data=darray, columns=names_out)
		
		df_out.to_csv('Output/grid_' + str(int(step)).zfill(4) + '.csv')

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

	def selectCell(self):

		i = random.randint(0,len(self.global_grid)-1)
		return(i)

	def turnover(self):

		i = self.selectCell()

		r = random.uniform(0,1)
		disp_pool = np.array([])

		if(r < Pdisp):
			disp_pool = self.getNeighbours(i)
			print('dispersal')
		else:
			disp_pool = self.global_grid[i].populations
			print('local')

		old = random.randint(0,len(self.global_grid[i].populations)-1)

		s = random.uniform(0,1)

		if(r < Pdisp and s < Pspec):
			self.global_grid[i].populations[old] = self.MaxSpec + 1
			self.MaxSpec += 1
			print('speciation')
		else:
			new = random.randint(0,len(disp_pool)-1)
			self.global_grid[i].populations[old] = disp_pool[new]



	
if(os.path.exists("Output")):
	shutil.rmtree("Output")
os.mkdir("Output")

g = Grid()
g.fillGrid(5)
g.printGrid(0)

t1 = time.time()
for i in range(0,50000001):
	print(i)
	g.turnover()

	if(i%50000==0):
		g.updateSpecies()
		g.printGrid(i/50000)
t2 = time.time()
print(t2-t1)




