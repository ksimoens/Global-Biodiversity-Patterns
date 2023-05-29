import numpy as np
import random

from parameters import*

class Local():

	def __init__(self,lon,lat,i,hab):
		self.lon = lon
		self.lat = lat
		self.index = i
		if(AreaHabitat):
			self.Nloc_i = int(hab*Nloc)
		else:
			self.Nloc_i = Nloc

		self.populations = np.zeros(self.Nloc_i)

	def fillLocal(self,Nspec):
			
		for i in range(0,len(self.populations)):
			self.populations[i] = random.randint(0,Nspec-1)



