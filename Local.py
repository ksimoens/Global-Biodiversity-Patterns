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
	def __init__(self,lon,lat,i,hab,temp=0):
		self.lon = lon # longitude
		self.lat = lat # latitude
		self.index = i  # index in global_grid of Grid object

		# if habitat area is used
		# AreaHabitat / Nloc defined in parameters.py
		if(AreaHabitat):
			# scale the local number of population with the habitat area proxy
			self.Nloc_i = int(hab*Nloc)
		else:
			self.Nloc_i = Nloc

		# if temperature-dependent rates are active
		if(TempTurnover or TempSpeciation):
			self.temp = temp # temperature

		# initialise container for the populations
		self.populations = []

	# fill the local grid cell with populations
	# Population defined in Population.py
	def fillLocal(self):

		for i in range(0,self.Nloc_i):
			self.populations.append( Population(self.index,i) )

