
import numpy as np 
from parameters import*

class Population():

	def __init__(self,i_glo,i_loc,tmin=0,tmax=0):

		self.glob_index = i_glo
		self.loc_index = i_loc
		self.species = -1
		if(TempNiches):
			self.T_min = tmin
			self.T_max = tmax