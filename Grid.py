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
		
		df_grid = pd.read_csv('grid_test_scaling.csv',index_col=0)
		k = 0
		for i in range(0,Nlon):
			for j in range(0,Nlat):
				self.global_grid.append( Local(df_grid.iloc[k,1],df_grid.iloc[k,0],k,df_grid.iloc[k,2]) )
				#self.global_grid.append( Local(i,j,k) )
				k += 1

		self.Nspec = 0
		self.species = np.zeros(self.Nspec)
		self.MaxSpec = 0
		#self.lat_max = df_grid['Y_COORD'].abs().max()

	def fillGrid(self):
		
		for loc in self.global_grid:
			loc.fillLocal()
		#self.Nspec = Nspec
		#self.species = np.arange(0,Nspec)
		#self.MaxSpec = np.max(self.species)

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


