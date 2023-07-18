# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Simulation of Mechanistic Model
#####################################
# LOCAL grid cell
#####################################


# ------------ IMPORT FROM OTHER FILES ---------------------

from parameters import*
from Population import*

# ----------------------------------------------------------


# ---------------- LOCAL CLASS -----------------------------

class Local():

	# define constructor
	# the Local object holds all information for the local grid cell
	def __init__(self,lon,lat,i,act,up,low,lef,rig,uplef,uprig,lowlef,lowrig,tem,hab):
		self.lon = lon # longitude
		self.lat = lat # latitude
		self.index = i # index in global_grid of Grid object
		self.active = act # is cell active?
		self.upper = up # is upper nghb active?
		self.lower = low # is lower nghb active?
		self.left = lef # is left nghb active? 
		self.right = rig # is right nghb active?
		self.upper_left = uplef # is upper left nghb active?
		self.upper_right = uprig # is upper right nghb active?
		self.lower_left = lowlef # is lower left nghb active?
		self.lower_right = lowrig # is lower right nghb active?
		self.temp = tem # temperature
		self.habit = hab # habitat area proxy

		# if habitat area is used
		# AreaHabitat / HabSlope / Nloc defined in parameters.py
		if(AreaHabitat):
			if(self.active == 1):
				# scale the local number of population with the habitat area proxy
				self.Nloc_i = int((HabSlope-1)*Nloc*self.habit - Nloc*(HabSlope-2))
			else:
				self.Nloc_i = float('nan')
		else:
			self.Nloc_i = Nloc

		# initialise container for the populations
		self.populations = []

	# fill the local grid cell with populations
	# Population defined in Population.py
	def fillLocal(self):

		for i in range(0,self.Nloc_i):
			self.populations.append( Population(self.index,i) )