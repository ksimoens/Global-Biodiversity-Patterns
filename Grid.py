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
		
		df_grid = pd.read_csv('simulation_grid.csv')
		k = 0
		for i in range(0,Nlon):
			for j in range(0,Nlat):
				self.global_grid.append( Local(df_grid.loc[k, 'x'],df_grid.loc[k, 'y'],k,df_grid.loc[k, 'active'],df_grid.loc[k, 'upper'],df_grid.loc[k, 'lower'],df_grid.loc[k, 'left'],df_grid.loc[k, 'right'],df_grid.loc[k, 'upper_left'],df_grid.loc[k, 'upper_right'],df_grid.loc[k, 'lower_left'],df_grid.loc[k, 'lower_right'],df_grid.loc[k, 'temp'],df_grid.loc[k, 'prec']) )
				k += 1

		self.Nspec = 0
		self.species = np.zeros(self.Nspec)
		self.MaxSpec = 0
		self.temp_min = df_grid.loc[:,'temp'].abs().max()

	def fillGrid(self):
		
		for loc in self.global_grid:
			if(loc.active == 1):
				loc.fillLocal()
		
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

		lon = np.array([-1,0,1,-1,0,1,-1,0,1])
		lat = np.array([1,1,1,0,0,0,-1,-1,-1])

		for i in range(0,len(lat_list)):
			lon_list[i] = lon[i]
			lat_list[i] = lat[i]

			for j in range(0,self.Nspec):
				spec_list[i][j] = np.count_nonzero(self.global_grid[i].populations == self.species[j])


		darray = np.concatenate((lon_list, lat_list)).reshape((-1, 2), order='F')
		darray = np.concatenate( (darray,spec_list) ,axis=1)
		names_out = ['x','y'] + ['spec_' + str(int(ID)) for ID in self.species] 		
		df_out = pd.DataFrame(data=darray, columns=names_out)
		
		df_out.to_csv('Output/grid_' + str(int(step)).zfill(4) + '.csv')

	def getNeighbours(self,index):

		neighbours = np.array([])

		up = self.global_grid[index].upper
		low = self.global_grid[index].lower
		lef = self.global_grid[index].left
		rig = self.global_grid[index].right
		up_lef = self.global_grid[index].upper_left
		up_rig = self.global_grid[index].upper_right
		low_lef = self.global_grid[index].lower_left
		low_rig = self.global_grid[index].lower_right

		# - - -
		# - o x
		# - - -
		if(rig == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index+1].populations))
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
			
		# - - x
		# - o -
		# - - -
		if(up_rig == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon+1].populations))

		# x - -
		# - o -
		# - - -
		if(up_lef == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon-1].populations))
	#	else: nothing
			# _ _ _
			# - o -
			# - - -

	
		# 6
		# - - -
		# - o -
		# - x -
		if(low == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon].populations))

		# - - -
		# - o -
		# - - x
		if(low_rig == 1):
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon+1].populations))

		# - - - 
		# - o -
		# x - -
		if(low_lef == 1):
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


