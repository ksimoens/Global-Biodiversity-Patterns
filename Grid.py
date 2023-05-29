import numpy as np
import pandas as pd
import random
import time
import os
import shutil

from parameters import*
from Local import*
from Species import*

class Grid():

	def __init__(self):

		self.global_grid = []
		
		df_grid = pd.read_csv('grid_test_scaling.csv',index_col=0)
		k = 0
		for i in range(0,Nlon):
			for j in range(0,Nlat):
				self.global_grid.append( Local(df_grid.iloc[k,1],df_grid.iloc[k,0],k,df_grid.iloc[k,2]) )
				k += 1

		if(TempSpeciation):
			lat_max = df_grid['Y_COORD'].abs().max()
			T_min = 303.15 - (1/3)*np.absolute(lat_max)
			self.Boltz_min = np.exp(-0.65 / 8.617e-5 / T_min)

		self.Nspec = 0
		self.species = []
		self.MaxSpec = 0
		self.probabilities = np.ones((len(self.global_grid)))

	def adjustProbabilities(self):
		for i in range(0,len(self.probabilities)):
			if(TempTurnover):
				T = 303.15 - (1/3)*np.absolute(self.global_grid[i].lat)
				self.probabilities[i] *= np.exp(-0.65 / 8.617e-5 / T)
			if(AreaHabitat):
				self.probabilities[i] *= self.global_grid[i].Nloc_i

		self.probabilities = self.probabilities / np.sum(self.probabilities,axis=0)

	def fillGrid(self,Nspec):
		
		if(TempNiches):
			for ilat in range(0,Nlat):
				 T = 303.15 - (1/3)*np.absolute(self.global_grid[ilat*Nlon].lat)
				 spec = Species(T,self.Nspec)
				 self.Nspec += 1
				 self.species.append(spec)
				 for ilon in range(0,Nlon):
				 	for i in range(0,len(self.global_grid[ilat*Nlon + ilon].populations)):
				 		self.global_grid[ilat*Nlon + ilon].populations[i] = spec.index

			self.MaxSpec = self.Nspec - 1

		else:
			for loc in self.global_grid:
				loc.fillLocal(Nspec)
			self.Nspec = Nspec
			self.species = np.arange(0,Nspec)
			self.MaxSpec = np.max(self.species)

	def updateSpecies(self):

		pop_list = np.array([])
		for cell in self.global_grid:
			pop_list = np.concatenate((pop_list,cell.populations))

		if(TempNiches):
			new_spec_list = []
			for spec in self.species:
				n = np.count_nonzero(pop_list == spec.index)
				if(n != 0):
					new_spec_list.append(spec)

			self.species = new_spec_list


		else:
			self.species = np.array([])
			for i in range(0,self.MaxSpec+1):
				n = np.count_nonzero(pop_list == i)
				if(n != 0):
					self.species = np.append(self.species,i)

		self.Nspec = len(self.species)
		print('Number of active species: ',str(self.Nspec))

	def printGrid(self,step):

		lon_list = np.zeros(len(self.global_grid))
		lat_list = np.zeros(len(self.global_grid))

		spec_list = np.zeros( (len(self.global_grid),self.Nspec) )

		for i in range(0,len(lat_list)):
			lon_list[i] = self.global_grid[i].lon
			lat_list[i] = self.global_grid[i].lat

			for j in range(0,self.Nspec):
				if(TempNiches):
					spec_list[i][j] = np.count_nonzero(self.global_grid[i].populations == self.species[j].index)
				else:
					spec_list[i][j] = np.count_nonzero(self.global_grid[i].populations == self.species[j])


		darray = np.concatenate((lon_list, lat_list)).reshape((-1, 2), order='F')
		darray = np.concatenate( (darray,spec_list) ,axis=1)
		if(TempNiches):
			ID_list_spec = np.array([])
			for spec in self.species:
				ID_list_spec = np.append(ID_list_spec,spec.index)
			names_out = ['lon','lat'] + ['spec_' + str(int(ID)) for ID in ID_list_spec]
		else: 
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

		if(TempTurnover):
			r = random.uniform(0,1)
			s = 0
			for i in range(0,len(self.probabilities)):
				s += self.probabilities[i]
				if(s > r):
					index = i
					break
		else:
			index = random.randint(0,len(self.global_grid)-1)
		
		return(index)

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

		old = random.randint(0,self.global_grid[i].Nloc_i-1)

		s = random.uniform(0,1)

		if(TempSpeciation):
			T = 303.15 - (1/3)*np.absolute(self.global_grid[i].lat)
			Pspec_i = Pspec * np.exp(-0.65 / 8.617e-5 / T) / self.Boltz_min
		else:
			Pspec_i = Pspec

		if(s < Pspec_i):
			self.global_grid[i].populations[old] = self.MaxSpec + 1
			if(TempNiches):
				T = 303.15 - (1/3)*np.absolute(self.global_grid[i].lat)
				self.species.append(Species(T,self.MaxSpec + 1))
			self.MaxSpec += 1
			print('speciation')
		else:
			valid = False
			while(not valid):
				new = random.randint(0,len(disp_pool)-1)
				if(TempNiches):
					T = 303.15 - (1/3)*np.absolute(self.global_grid[i].lat)
					i_spec = disp_pool[new]
					new_spec = None
					for spec in self.species:
						if(spec.index == i_spec):
							new_spec = spec
							break
					if(T > new_spec.min_temp and T < new_spec.max_temp):
						valid = True

				else:
					valid = True

			self.global_grid[i].populations[old] = disp_pool[new]



	
if(os.path.exists("Output")):
	shutil.rmtree("Output")
os.mkdir("Output")

g = Grid()
g.fillGrid(1)
g.adjustProbabilities()

t1 = time.time()
for i in range(0,5000001):
	print(i)
	g.turnover()

	if(i%10000==0):
		g.updateSpecies()
		g.printGrid(i/10000)
t2 = time.time()
print(t2-t1)



