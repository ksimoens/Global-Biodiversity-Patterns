# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Simulation of Mechanistic Model
#####################################
# local POPULATION
#####################################


# ----------------- IMPORT MODULES -------------------------

import numpy as np

# ----------------------------------------------------------


# ---------------- POPULATION CLASS ------------------------

class Population():

	# define constructor
	def __init__(self,i_glo,i_loc):

		self.glob_index = i_glo # index in global_grid (Grid)
		self.loc_index = i_loc # index in the local cell (Local)
		self.species = -1 # species identifier; unidentified = -1