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
import multiprocessing as mp

# ----------------------------------------------------------


# ---------------- SPECIES CLASS ---------------------------

class Species():

	# define constructor
	def __init__(self,o,i_glo,i_loc,T_min=-1,T_max=-1):
		# unique species identifier = order in which it appeared in the simulation
		self.order = o
		# global index of the cell in which the species appeared
		self.glob_index = i_glo
		# local index of the population in which the species appeared
		self.loc_index = i_loc
		# if temperature ranges are active (defined in parameters.py)
		if(TempNiches):
			# lower and upper end of the range in which the central temperature can lie
			self.temp_min = T_min
			self.temp_max = T_max

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
# very ugly; has to be rewritten and compartmentalised

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
	# if temperature ranges are active (defined in parameters.py)
	if(TempNiches):
		# containers for limits of the range for the central temperature
		temp_min_list = []
		temp_max_list = []
	# container for the temperature-dependent rates (defined in parameters.py)
	# not normalised for speciation rate scaling
	prob_list = []
	for i in range(0,len(g.global_grid)):
		for j in range(0,len(g.global_grid[i].populations)):
			if(g.global_grid[i].active == 1):
				glob_index.append(i)
				loc_index.append(j)
				if(TempTurnover):	
					T = g.global_grid[i].temp
					# Boltzmann factor
					prob_list.append(np.exp(-0.65 / 8.617e-5 / T))
				else:
					prob_list.append(-1)
				if(TempNiches):
					T = g.global_grid[i].temp
					temp_min_list.append(T)
					temp_max_list.append(T)

	if(TempSpeciation):
		# Boltzmann factor for the lowest temperature
		boltz_min = np.exp(-0.65 / 8.617e-5 / g.temp_min)

	# create the coalescence tree = dataframe
	# rows = individual populations
	if(TempNiches):
		tree = pd.DataFrame(data=np.vstack( (glob_index,loc_index,prob_list,temp_min_list,temp_max_list)).T,columns=['glob','loc','prob','temp_min','temp_max'])
	else:
		tree = pd.DataFrame(data=np.vstack( (glob_index,loc_index,prob_list)).T,columns=['glob','loc','prob'])
	tree = tree.astype({'glob':'int64'})
	tree = tree.astype({'loc':'int64'})

	# create reference IDlist for bookkeeping
	# rows = individual populations
	if(TempNiches):
		IDlist = pd.DataFrame(data=np.vstack( (glob_index,loc_index,temp_min_list,temp_max_list)).T,columns=['glob','loc','temp_min','temp_max'])
	else:
		IDlist = pd.DataFrame(data=np.vstack( (glob_index,loc_index)).T,columns=['glob','loc'])
	# container column for the species identity of the population
	IDlist['species'] = -1

	# counter for the steps in the coalescence algorithm
	count = 0
	# initial size of the coalescence tree
	tree_0 = len(str(len(tree)))

	# as long as more than one population remains, continue the algorithm
	while(len(tree) > 1):

		print("iteration:" + '\t' + str(count).zfill(8) + '\t' + 'tree size:' + '\t' + str(len(tree)).zfill(tree_0), end='\r' )
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
			# migration: pool of dispersers = Moore neighbourhood
			disp_pool = g.getNeighbours(old_pop.glob_index)
		else:
			# no migration: pool of dispersors = local populations
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

		# if temperature ranges are active
		if(TempNiches):
			# is migration allowed?
			valid = False
			# count how many migration attempts happened
			vcount = 0
			# as long as the migration is not allowed, repeat the algorithm
			while(not valid):

				# randomly select a population from the disperser pool
				rNew = random.randint(0,len(disp_pool)-1)
				new_pop = disp_pool[rNew]

				# get environmental temperature values at both locations
				T_old = g.global_grid[old_pop.glob_index].temp
				T_new = g.global_grid[new_pop.glob_index].temp

				# range of the central temperature of the common ancestor
				temp_ancest_min = max(old_pop.T_min,new_pop.T_min) 
				temp_ancest_max = min(old_pop.T_max,new_pop.T_max)

				# if speciation
				if(rSpec < Pspec_i):

					# if dispersal is allowed
					if(temp_ancest_max > temp_ancest_min and T_old < temp_ancest_max and T_old > temp_ancest_min and T_new < temp_ancest_max and T_new > temp_ancest_min):
						
						# the new temperature limits of the disperser = temperature limits of the common ancestor
						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_min = temp_ancest_min
						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_max = temp_ancest_max
						# the new temperature limits of the deceased = default range
						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_min = T_old - NicheWidth
						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_max = T_old + NicheWidth
						# dispersal is allowed
						valid = True
						# create the new species and add to the list
						new_spec = Species(len(spec_list),old_pop.glob_index,old_pop.loc_index,T_old-NicheWidth,T_old+NicheWidth)
						spec_list.append(new_spec)
						# identify deceased in the ID list
						IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index) & (IDlist['species']==-1), 'species'] = new_spec.order
						# remove the deceased from the tree
						tree = tree.drop(tree.index[r])

				# if no speciation, but coalescence
				elif(len(tree[(tree['glob']==new_pop.glob_index) & (tree['loc']==new_pop.loc_index)]) != 0):
					
					# if dispersal is allowed
					if(temp_ancest_max > temp_ancest_min and T_old < temp_ancest_max and T_old > temp_ancest_min and T_new < temp_ancest_max and T_new > temp_ancest_min):

						# the new temperature limits of the deceased = default range
						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_min = T_old - NicheWidth
						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_max = T_old + NicheWidth
						# the new temperature limits of the disperser = temperature limits of the common ancestor
						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_min = temp_ancest_min
						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_max = temp_ancest_max
						# get a copy from the index ID lists
						glob_list = copy.deepcopy(IDlist['glob'])
						loc_list = copy.deepcopy(IDlist['loc'])
						# replace the indices of the deceased with the ones from the disperser in the ID list
						IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'glob' ] = new_pop.glob_index
						IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'loc' ] = new_pop.loc_index
						# remove the deceased from the tree
						tree = tree.drop(tree.index[r])
						# dispersal is allowed
						valid = True

				# if no speciation and no coalescence
				else:
					
					# if dispersal is allowed		
					if(temp_ancest_max > temp_ancest_min and T_old < temp_ancest_max and T_old > temp_ancest_min and T_new < temp_ancest_max and T_new > temp_ancest_min):

						# the new temperature limits of the deceased = default range
						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_min = T_old - NicheWidth
						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_max = T_old + NicheWidth
						# the new temperature limits of the disperser = temperature limits of the common ancestor
						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_min = temp_ancest_min
						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_max = temp_ancest_max
						# dispersal is allowed
						valid = True
						# get a copy from the index ID lists
						glob_list = copy.deepcopy(IDlist['glob'])
						loc_list = copy.deepcopy(IDlist['loc'])
						# replace the indices of the deceased with the ones from the disperser in the ID list and the coalescence tree
						IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'glob' ] = new_pop.glob_index
						IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'loc' ] = new_pop.loc_index
						tree.iloc[r,tree.columns.get_loc('glob')] = new_pop.glob_index
						tree.iloc[r,tree.columns.get_loc('loc')] = new_pop.loc_index
						# if temperature-dependent rates are active (defined in parameters.py)
						if(TempTurnover):
							# adjust the probability to the new temperature
							tree.iloc[r,tree.columns.get_loc('prob')] = np.exp(-0.65 / 8.617e-5 / g.global_grid[new_pop.glob_index].temp)
				
				# if dispersal should not be allowed but it is by the algorithm: stop the simulation
				if( (new_pop.T_min > T_old or new_pop.T_max < T_old) and valid == True):
					print('T_old ' + str(g.global_grid[old_pop.glob_index].temp))
					print('T_old_min ' + str(old_pop.T_min))
					print('T_old_max ' + str(old_pop.T_max))
					print('T_new_min ' + str(new_pop.T_min))
					print('T_new_max ' + str(new_pop.T_max))
					input()
				
				# number of dispersal trials increases
				vcount += 1

		# if temperature ranges are inactive
		else:

			# randomly select a population from the disperser pool
			rNew = random.randint(0,len(disp_pool)-1)
			new_pop = disp_pool[rNew]

			# if speciation 
			if(rSpec < Pspec_i):
				# create new species and add to the list
				new_spec = Species(len(spec_list),old_pop.glob_index,old_pop.loc_index)
				spec_list.append(new_spec)
				# identify deceased in the ID list
				IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index) & (IDlist['species']==-1), 'species'] = new_spec.order
				# remove the deceased from the list
				tree = tree.drop(tree.index[r])

			# if no speciation
			else:
				# if coalescence
				if len(tree[(tree['glob']==new_pop.glob_index) & (tree['loc']==new_pop.loc_index)]) != 0:
					# remove the deceased from the tree
					tree = tree.drop(tree.index[r])
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
		if(not SpeciesInformation):
			# container for species identifier for each local population
			spec_sublist = np.zeros(len(g.global_grid[i].populations))
		for j in range(0,len(g.global_grid[i].populations)):
			# identify populations in the grid
			g.global_grid[i].populations[j].species = IDlist['species'][c]
			if(not SpeciesInformation):
				# identify populations in the container
				spec_sublist[j] = IDlist['species'][c]
			c += 1
		if(not SpeciesInformation):
			# calculate the diversity of the grid cell
			diver_array[i] = len(np.unique(spec_sublist))

	# create containers for cell coordinates
	lon_list = np.zeros(len(g.global_grid))
	lat_list = np.zeros(len(g.global_grid))

	# create container for the diversity matrix
	if(SpeciesInformation):
		spec_array = np.zeros( (len(g.global_grid),len(spec_list)) )

	# for each grid cell
	for i in range(0,len(lat_list)):
		# get coordinates
		lon_list[i] = g.global_grid[i].lon
		lat_list[i] = g.global_grid[i].lat

		if(SpeciesInformation):
			# for each species in the species list
			for j in range(0,len(spec_list)):
				n = 0
				# count how many populations belong to the species
				for x in range(0,len(g.global_grid[i].populations)):
					if(g.global_grid[i].populations[x].species == spec_list[j].order):
						n += 1	
				spec_array[i][j] = n

	# if the diversity matrix is active (defined in parameters.py)
	if(SpeciesInformation):
		# create species names (= identifiers)
		spec_ID = np.zeros(len(spec_list))
		for i in range(0,len(spec_ID)):
			spec_ID[i] = spec_list[i].order

		# create output dataframe with diversity matrix 
		darray = np.concatenate((lon_list, lat_list)).reshape((-1, 2), order='F')
		darray = np.concatenate( (darray,spec_array) ,axis=1)
		names_out = ['lon','lat'] + ['spec_' + str(int(ID)) for ID in spec_ID] 		
		df_out = pd.DataFrame(data=darray, columns=names_out)
	# if not active
	else:
		# create output dataframe with diversity column
		darray = np.concatenate((lon_list,lat_list,diver_array)).reshape((-1, 3), order='F')
		names_out = ['lon','lat','diversity']
		df_out = pd.DataFrame(data=darray, columns=names_out)

	# write out the output
	# identifier = replicate number
	df_out.to_csv('Output/grid_' + str(int(k)).zfill(4) + '.csv')