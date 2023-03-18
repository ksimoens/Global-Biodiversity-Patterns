import numpy as np
import random

from parameters import*
from Population import*

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
		self.populations = []

	def fillLocal(self):

		for i in range(0,self.Nloc_i):
			if(TempNiches):
				T = 303.15 - (1/3)*np.absolute(self.lat)
				self.populations.append( Population(self.index,i,T - NicheWidth, T + NicheWidth) )
			else:
				self.populations.append( Population(self.index,i) )

