# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Simulation of Mechanistic Model
#####################################
# coalescence TREE
#####################################


# ------------ IMPORT FROM OTHER FILES ---------------------

from parameters import*
from Grid import*

# ----------------------------------------------------------


# ----------------- IMPORT MODULES -------------------------

import numpy as np
import pandas as pd
import random
import copy
import time

# ----------------------------------------------------------


# ---------------- SPECIES CLASS ---------------------------

class Species():

	# define constructor
	def __init__(self,o,i_glo,i_loc):
		# unique species identifier = order in which it appeared in the simulation
		self.order = o
		# global index of the cell in which the species appeared
		self.glob_index = i_glo
		# local index of the population in which the species appeared
		self.loc_index = i_loc

# ----------------------------------------------------------


# ----------------- SELECT FOR DEATH -----------------------
# select a cell to be disturbed according to the rules defined in parameters.py
# requires a coalescence tree as input

def selectCell(t):
	# if temperature-dependent death rates are active
	if(TempTurnover):
		r = random.uniform(0,1)
		s = 0
		index = -1
		# select a cell based on a temperature-scaled probability 
		# normalise
		probs = t['prob'].to_numpy()/np.sum(t['prob'])
		for i in range(0,len(t)):
			s += probs[i]
			if(s > r):
				index = i
				break
	# if temperature-dependent death rates are inactive
	else:
		index = random.randint(0,len(t)-1)

	return(int(index))

# ----------------------------------------------------------


# --------------- COALESCENCE SIMULATION -------------------
# the working horse of the simulation
# runs the entire coalescence algorithm
# very ugly; has to be rewritten and compartimentalised

# requires a replicate identifier (k) as input 

def replicateRun(k):

	print("---------------------------------------------")
	print("SIMULATION RUN " + str(k+1) + " of " + str(Nrep))
	print("---------------------------------------------")

	# create and initialise grid
	g = Grid()
	g.fillGrid()

	# create the initial coalescence tree and population identifier

	# species identifier container
	spec_list = []
	# container for the index in the global_grid (Grid)
	glob_index = []
	# container for the index in the local cell (Local)
	loc_index = []
	# container for the temperature-dependent rates (defined in parameters.py)
	# not normalised for speciation rate scaling
	prob_list = []
	for i in range(0,len(g.global_grid)):
		for j in range(0,len(g.global_grid[i].populations)):
			glob_index.append(i)
			loc_index.append(j)
			if(TempTurnover):	
				T = g.global_grid[i].temp
				prob_list.append(np.exp(-0.65 / 8.617e-5 / T))
			else:
				prob_list.append(-1)


	if(TempSpeciation):
		# Boltzmann factor for the lowest temperature
		boltz_min = np.exp(-0.65 / 8.617e-5 / g.Tmin)

	# create the coalescence tree = dataframe
	# rows = individual populations
	tree = pd.DataFrame(data=np.vstack( (glob_index,loc_index,prob_list)).T,columns=['glob','loc','prob'])
	tree = tree.astype({'glob':'int64'})
	tree = tree.astype({'loc':'int64'})

	# create reference IDlist for bookkeeping
	# rows = individual populations
	IDlist = pd.DataFrame(data=np.vstack( (glob_index,loc_index)).T,columns=['glob','loc'])
	# container column for the species identity of the population
	IDlist['species'] = -1

	# counter for the steps in the coalescence algorithm
	count = 0
	# count number of coalescence / speciation / migration events
	count_coal = 0
	count_spec = 0
	count_migr = 0

	# set seed of random number generator for all child processes
	pid = os.getpid()
	random.seed(pid*time.time())
	# initial size of the coalescence tree
	tree_0 = len(str(len(tree)))

	# as long as more than one population remains, continue the algorithm
	while(len(tree) > 1):

		print("iteration:" + '\t' + str(count).zfill(10) + '\t' + 'tree size:' + '\t' + str(len(tree)).zfill(tree_0), end='\r' )
		# select an index for a population in the tree to die
		r = selectCell(tree)
		# population selected to die in grid
		old_pop = g.global_grid[tree['glob'].iloc[r]].populations[tree['loc'].iloc[r]]

		# random number for migration
		rDisp = random.uniform(0,1)
		# initialise pool of possible dispersers
		disp_pool = []
		# check for migration (Pdisp defined in parameters.py)
		if rDisp < Pdisp:
			# alternative migration (defined in parameters.py)
			if(NewMigration):
				disp_pool = g.getAllPopulations(old_pop.glob_index)
			else:
				disp_pool = g.getNeighbours(old_pop.glob_index)
			count_migr += 1
		else:
			if(NewMigration):
				disp_pool = g.getNeighbours(old_pop.glob_index)
				disp_pool = np.concatenate((disp_pool,g.global_grid[old_pop.glob_index].populations))
				disp_pool = np.delete(disp_pool,np.where(disp_pool == old_pop))
			else:
				disp_pool = g.global_grid[old_pop.glob_index].populations
				disp_pool = np.delete(disp_pool,old_pop.loc_index)
		
		# get the speciation probability for the cell
		# Pspec / TempSpeciation defined in parameters.py
		if(TempSpeciation):
			Pspec_i = Pspec * tree['prob'].iloc[r] / boltz_min
		else:
			Pspec_i = Pspec

		# random number for speciation
		rSpec = random.uniform(0,1)

		# if speciation 
		if(rSpec < Pspec_i):
			count_spec += 1
			count_coal += 1
			# create new species and add to the list
			new_spec = Species(len(spec_list),old_pop.glob_index,old_pop.loc_index)
			spec_list.append(new_spec)
			# identify deceased in the ID list
			IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index) & (IDlist['species']==-1), 'species'] = new_spec.order
			# remove the deceased from the list
			tree = tree.drop(tree.index[r])

		# if no speciation
		else:
			# randomly select a population from the disperser pool
			rNew = random.randint(0,len(disp_pool)-1)
			new_pop = disp_pool[rNew]
			
			# if coalescence
			if len(tree[(tree['glob']==new_pop.glob_index) & (tree['loc']==new_pop.loc_index)]) != 0:
				# remove the deceased from the tree
				tree = tree.drop(tree.index[r])
				count_coal += 1
			# if no coalescence
			else:
				# replace the indices of the deceased with the ones from the disperser in the coalescence tree
				tree.iloc[r,tree.columns.get_loc('glob')] = new_pop.glob_index
				tree.iloc[r,tree.columns.get_loc('loc')] = new_pop.loc_index
				# if temperature-dependent rates are active (defined in parameters.py)
				if(TempTurnover):
					# adjust the probability to the new temperature
					T = g.global_grid[new_pop.glob_index].temp
					tree.iloc[r,tree.columns.get_loc('prob')] = np.exp(-0.65 / 8.617e-5 / T)

			# get a copy from the index ID lists
			glob_list = copy.deepcopy(IDlist['glob'])
			loc_list = copy.deepcopy(IDlist['loc'])
			# replace the indices of the deceased with the ones from the disperser in the ID list
			IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'glob' ] = new_pop.glob_index
			IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'loc' ] = new_pop.loc_index

		# number of coalescence steps increases
		count += 1

	# handle the final population in the tree
	final_pop = g.global_grid[tree['glob'].iloc[0]].populations[tree['loc'].iloc[0]]
	# create a new species for the final population and add to the list
	final_spec = Species(len(spec_list),tree.iloc[0]['glob'],tree.iloc[0]['loc'])
	spec_list.append(final_spec)
	# identify final population in the ID list
	IDlist.loc[ (IDlist['glob']==tree['glob'].iloc[0]) & (IDlist['loc']==tree['loc'].iloc[0]) & (IDlist['species']==-1), 'species'] = final_spec.order

	# create the diversity column if the diversity matrix is not active (defined in parameters.py)
	if(not SpeciesInformation):
		diver_array = np.zeros(len(g.global_grid))

		# index in the ID list
		c = 0
		# go over all cells in the grid
		for i in range(0,len(g.global_grid)):
			# container for species identifier for each local population
			spec_sublist = np.zeros(len(g.global_grid[i].populations))
			for j in range(0,len(g.global_grid[i].populations)):
				# identify populations in the grid
				spec_sublist[j] = IDlist['species'][c]
				c += 1

			# calculate the diversity of the grid cell
			diver_array[i] = len(np.unique(spec_sublist))

		# create containers for cell coordinates
		lon_list = np.zeros(len(g.global_grid))
		lat_list = np.zeros(len(g.global_grid))

		# for each grid cell
		for i in range(0,len(lat_list)):
			# get coordinates
			lon_list[i] = g.global_grid[i].lon
			lat_list[i] = g.global_grid[i].lat

	# if species information is active (defined in parameters.py)
	if(SpeciesInformation):	
		# get a list of species identities for each population
		# a lot faster than writing out diversity matrices
		df_out = pd.DataFrame(data=IDlist['species'],columns=['species'])
	# if not active
	else:
		# create output dataframe with diversity column
		darray = np.concatenate((lon_list,lat_list,diver_array)).reshape((-1, 3), order='F')
		names_out = ['lon','lat','diversity']
		df_out = pd.DataFrame(data=darray, columns=names_out)

	# write out the output
	# identifier = replicate number	
	df_out.to_csv('Output/grid_' + str(int(k)).zfill(4) + '.csv')
	# some additional information on the simulation behaviour
	f = open("Output/count_" + str(int(k)).zfill(4) + ".txt", "x")
	f.write('count' + '\t' + str(count) + '\n')
	f.write('count_spec' + '\t' + str(count_spec) + '\n')
	f.write('count_coal' + '\t' + str(count_coal) + '\n')
	f.write('count_migr' + '\t' + str(count_migr) + '\n')
	f.close()