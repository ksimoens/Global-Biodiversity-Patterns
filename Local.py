import numpy as np
import random

from parameters import*

class Local():

	def __init__(self,lon,lat,i):
		self.lon = lon
		self.lat = lat
		self.index = i
		if(AreaHabitat):
			factor = -4./90.*np.absolute(lat) + 5.
			self.Nloc_i = int(factor*Nloc)
		else:
			self.Nloc_i = Nloc
		self.populations = np.zeros(self.Nloc_i)

	def fillLocal(self,Nspec):

		for i in range(0,len(self.populations)):
			self.populations[i] = random.randint(0,Nspec-1)

