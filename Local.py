import numpy as np
import random

from parameters import*
from Population import*

class Local():

	def __init__(self,lon,lat,i,act,up,low,lef,rig,uplef,uprig,lowlef,lowrig,tem,hab):
		self.lon = lon
		self.lat = lat
		self.index = i
		self.active = act
		self.upper = up
		self.lower = low
		self.left = lef
		self.right = rig
		self.upper_left = uplef
		self.upper_right = uprig
		self.lower_left = lowlef
		self.lower_right = lowrig
		self.temp = tem
		self.habit = hab
		if(AreaHabitat):
			if(self.active == 1):
				self.Nloc_i = int((HabSlope-1)*Nloc*self.habit - Nloc*(HabSlope-2))
			else:
				self.Nloc_i = float('nan')
		else:
			self.Nloc_i = Nloc
		self.populations = []

	def fillLocal(self):

		for i in range(0,self.Nloc_i):
			self.populations.append( Population(self.index,i) )