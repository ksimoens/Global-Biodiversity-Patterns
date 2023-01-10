import numpy as np
import pandas as pd

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

	def fillGrid(self,Nspec):
		
		for loc in self.global_grid:
			loc.fillLocal(Nspec)
		self.Nspec = Nspec
		self.species = np.arange(0,Nspec)

	def printGrid(self):

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
		names_out = ['lon','lat'] + ['spec_' + str(ID) for ID in self.species] 		
		df_out = pd.DataFrame(data=darray, columns=names_out)
		
		df_out.to_csv('grid.csv')

	def getNeighbours(self,index):

		neighbours = np.array([])

		# 1
		# | - - . - - |
		# | x - . - o |
		# | - - . - - |
		if((index+1) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon+1].populations))
			print(index-Nlon+1)
		# - - -
		# - o x
		# - - -
		else:
			neighbours = np.concatenate((neighbours , self.global_grid[index+1].populations))
			print(index+1)
		# 2
		# | - - . - - |
		# | o - . - x |
		# | - - . - - |
		if((index) % Nlon == 0):
			neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon-1].populations))
			print(index+Nlon-1)
		# - - -
		# x o -
		# - - -
		else:
			neighbours = np.concatenate((neighbours , self.global_grid[index-1].populations))
			print(index-1)

		if(index >= Nlon):
			# 3
			# - x -
			# - o -
			# - - -
			neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon].populations))
			print(index-Nlon)
			
			# 4
			# | x - . - - |
			# | - - . - o |
			# | - - . - - |
			if((index+1) % Nlon == 0):
				neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon-Nlon+1].populations))
				print(index-Nlon-Nlon+1)
			# - - x
			# - o -
			# - - -
			else:
				neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon+1].populations))
				print(index-Nlon+1)

			# 5
			# | - - . - x |
			# | o - . - - |
			# | - - . - - |
			if((index) % 45 == 0):
				neighbours = np.concatenate((neighbours , self.global_grid[index-1].populations))
				print(index-1)
			# x - -
			# - o -
			# - - -
			else:
				neighbours = np.concatenate((neighbours , self.global_grid[index-Nlon-1].populations))
				print(index-Nlon-1)
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
			print(index+Nlon)

			# 7
			# | - - . - - |
			# | - - . - o |
			# | x - . - - |
			if((index+1) % Nlon == 0):
				neighbours = np.concatenate((neighbours , self.global_grid[index+1].populations))
				print(index+1)
			# - - -
			# - o -
			# - - x
			else:
				neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon+1].populations))
				print(index+Nlon+1)

			# 8
			# | - - . - - |
			# | o - . - - |
			# | - - . - x |
			if((index) % 45 == 0):
				neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon+Nlon-1].populations))
				print(index+Nlon+Nlon-1)
			# - - - 
			# - o -
			# x - -
			else:
				neighbours = np.concatenate((neighbours , self.global_grid[index+Nlon-1].populations))
				print(index+Nlon-1)
			
			
	#	else: nothing
			# - - -
			# - o -
			# _ _ _



	


g = Grid()
g.fillGrid(5)
#g.printGrid()
g.getNeighbours(584)

