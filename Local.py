import numpy as np
import random

from parameters import*

class Local():

	def __init__(self,lon,lat,i):
		self.populations = np.zeros(Nloc)
		self.lon = lon
		self.lat = lat
		self.index = i

	def fillLocal(self,Nspec):

		for i in range(0,len(self.populations)):
			self.populations[i] = int(random.randint(0,Nspec))

