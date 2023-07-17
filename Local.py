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

# ----------------------------------------------------------


# ----------------- IMPORT MODULES -------------------------

import numpy as np
import random

# ----------------------------------------------------------


# ---------------- LOCAL CLASS -----------------------------
class Local():

	# define constructor
	def __init__(self,lon,lat,i):
		# initialise container for the populations
		# number of local populations defined in parameters.py
		self.populations = np.zeros(Nloc)
		# set the grid cell coordinates
		self.lon = lon
		self.lat = lat
		# set the index in the global grid
		self.index = i

	# fill the local grid cell with populations
	# requires an initial number of species: Nspec
	def fillLocal(self,Nspec):

		# fill with integers = species identifiers
		for i in range(0,len(self.populations)):
			self.populations[i] = random.randint(0,Nspec-1)

