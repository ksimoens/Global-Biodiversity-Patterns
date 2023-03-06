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

		self.global_grid = np.empty([int(Nlat*np.sqrt(Nloc)),int(Nlon*np.sqrt(Nloc))])
		
		
		#df_grid = pd.read_csv('land_all_data_inc_working.csv').iloc[:,[0,5,6]]

		self.Nspec = Nspec0
		self.species = np.zeros(self.Nspec)
		self.MaxSpec = Nspec0-1 
		print(self.global_grid.shape[0])

	def fillGrid(self):
		
		c = 0
		for i in range(0,Nlat):
			for j in range(0,Nlon):
				self.global_grid[i][j] = c#random.randint(0,Nspec0-1)
				c += 1



	def updateSpecies(self):

		pop_list = np.array([])
		for pop in self.global_grid:
			pop_list = np.append((pop_list,pop))

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

		lon = np.array([-1,0,1,-1,0,1,-1,0,1])
		lat = np.array([1,1,1,0,0,0,-1,-1,-1])

		for i in range(0,len(lat_list)):
			lon_list[i] = lon[i]
			lat_list[i] = lat[i]

			for j in range(0,self.Nspec):
				spec_list[i][j] = np.count_nonzero(self.global_grid[i].populations == self.species[j])


		darray = np.concatenate((lon_list, lat_list)).reshape((-1, 2), order='F')
		darray = np.concatenate( (darray,spec_list) ,axis=1)
		names_out = ['lon','lat'] + ['spec_' + str(int(ID)) for ID in self.species] 		
		df_out = pd.DataFrame(data=darray, columns=names_out)
		
		df_out.to_csv('Output/grid_' + str(int(step)).zfill(4) + '.csv')

	def getNeighbours(self,i_lat,i_lon):

		neighbours = np.array([])

		# 1
		# | - - . - - |
		# | x - . - o |
		# | - - . - - |
		if(i_lon == Nlon-1):
			neighbours = np.append(neighbours , self.global_grid[i_lat][0])
		# - - -
		# - o x
		# - - -
		else:
			neighbours = np.append(neighbours , self.global_grid[i_lat][i_lon+1])
		# 2
		# | - - . - - |
		# | o - . - x |
		# | - - . - - |
		if(i_lon == 0):
			neighbours = np.append(neighbours , self.global_grid[i_lat][Nlon-1])
		# - - -
		# x o -
		# - - -
		else:
			neighbours = np.append(neighbours , self.global_grid[i_lat][i_lon-1])

		if(i_lat > 0):
			# 3
			# - x -
			# - o -
			# - - -
			neighbours = np.append(neighbours , self.global_grid[i_lat-1][i_lon])
			
			# 4
			# | x - . - - |
			# | - - . - o |
			# | - - . - - |
			if(i_lon == Nlon-1):
				neighbours = np.append(neighbours , self.global_grid[i_lat-1][0])
			# - - x
			# - o -
			# - - -
			else:
				neighbours = np.append(neighbours , self.global_grid[i_lat-1][i_lon+1])

			# 5
			# | - - . - x |
			# | o - . - - |
			# | - - . - - |
			if(i_lon == 0):
				neighbours = np.append(neighbours , self.global_grid[i_lat-1][Nlon-1])
			# x - -
			# - o -
			# - - -
			else:
				neighbours = np.append(neighbours , self.global_grid[i_lat-1][i_lon-1])
	#	else: nothing
			# _ _ _
			# - o -
			# - - -

		if(i_lat < Nlat-1):
			# 6
			# - - -
			# - o -
			# - x -
			neighbours = np.append(neighbours , self.global_grid[i_lat+1][i_lon])

			# 7
			# | - - . - - |
			# | - - . - o |
			# | x - . - - |
			if(i_lon == Nlon-1):
				neighbours = np.append(neighbours , self.global_grid[i_lat+1][0])
			# - - -
			# - o -
			# - - x
			else:
				neighbours = np.append(neighbours , self.global_grid[i_lat+1][i_lon+1])

			# 8
			# | - - . - - |
			# | o - . - - |
			# | - - . - x |
			if(i_lon == 0):
				neighbours = np.append(neighbours , self.global_grid[i_lat+1][Nlon-1])
			# - - - 
			# - o -
			# x - -
			else:
				neighbours = np.append(neighbours , self.global_grid[i_lat+1][i_lon-1])
			
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

		if(s < Pspec):
			self.global_grid[i].populations[old] = self.MaxSpec + 1
			self.MaxSpec += 1
			print('speciation')
		else:
			new = random.randint(0,len(disp_pool)-1)
			self.global_grid[i].populations[old] = disp_pool[new]


'''
	
if(os.path.exists("Output")):
	shutil.rmtree("Output")
os.mkdir("Output")

g = Grid()
g.fillGrid(1)
g.printGrid(0)


t1 = time.time()
for i in range(0,100001):
	print(i)
	g.turnover()

	if(i%100==0):
		g.updateSpecies()
		g.printGrid(i/100)
t2 = time.time()
print(t2-t1)

'''


