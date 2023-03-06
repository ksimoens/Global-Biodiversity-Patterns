import numpy as np
import random

from parameters import*
from Population import*

class Local():

	def __init__(self,lon,lat):
		self.populations = []
		self.lon = lon
		self.lat = lat

	def fillLocal(self,N):

		for i in range(0,N):
			self.populations.append(Population(self.index,i))

