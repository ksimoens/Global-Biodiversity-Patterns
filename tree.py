from parameters import*
from Grid import*

import numpy as np
import pandas as pd
import random
import copy
import time

class Species():

	def __init__(self,o,i_glo,i_loc,T_min=-1,T_max=-1):
		self.order = o
		self.glob_index = i_glo
		self.loc_index = i_loc
		if(TempNiches):
			self.temp_min = T_min
			self.temp_max = T_max

def selectCell(t):
	if(TempTurnover):
		r = random.uniform(0,1)
		s = 0
		index = -1
		probs = t['prob'].to_numpy()/np.sum(t['prob'])
		for i in range(0,len(t)):
			s += probs[i]
			if(s > r):
				index = i
				print(i)
				break
	else:
		index = random.randint(0,len(t)-1)

	return(int(index))

if(os.path.exists("Output")):
	shutil.rmtree("Output")
os.mkdir("Output")

t1 = time.time()

Nrange = 20

for k in range(0,Nrange):

	print("---------------------------------------------")
	print("SIMULATION RUN " + str(k) + " of " + str(Nrange))
	print("---------------------------------------------")

	g = Grid()
	g.fillGrid()

	spec_list = []

	glob_index = []
	loc_index = []

	if(TempNiches):
		temp_min_list = []
		temp_max_list = []

	prob_list = []
	for i in range(0,len(g.global_grid)):
		for j in range(0,len(g.global_grid[i].populations)):
			glob_index.append(i)
			loc_index.append(j)
			if(TempTurnover):	
				T = 303.15 - (1/3)*np.absolute(g.global_grid[i].lat)
				prob_list.append(np.exp(-0.65 / 8.617e-5 / T))
			else:
				prob_list.append(-1)
			if(TempNiches):
				T = 303.15 - (1/3)*np.absolute(g.global_grid[i].lat)
				temp_min_list.append(T)
				temp_max_list.append(T)

	#prob_list = prob_list / np.sum(prob_list)

	if(TempSpeciation):
		temp_min = 303.15 - (1/3)*np.absolute(g.lat_max)
		boltz_min = np.exp(-0.65 / 8.617e-5 / temp_min)

	if(TempNiches):
		tree = pd.DataFrame(data=np.vstack( (glob_index,loc_index,prob_list,temp_min_list,temp_max_list)).T,columns=['glob','loc','prob','temp_min','temp_max'])
	else:
		tree = pd.DataFrame(data=np.vstack( (glob_index,loc_index,prob_list)).T,columns=['glob','loc','prob'])
	tree = tree.astype({'glob':'int64'})
	tree = tree.astype({'loc':'int64'})

	if(TempNiches):
		IDlist = pd.DataFrame(data=np.vstack( (glob_index,loc_index,temp_min_list,temp_max_list)).T,columns=['glob','loc','temp_min','temp_max'])
	else:
		IDlist = pd.DataFrame(data=np.vstack( (glob_index,loc_index)).T,columns=['glob','loc'])
	IDlist['species'] = -1
	#print(IDlist)

	count = 0





	while(len(tree) > 1):

		print("iteration:" + '\t' + str(count) + '\t' + str(len(tree)), end='\r')
		r = selectCell(tree)
		old_pop = g.global_grid[tree['glob'].iloc[r]].populations[tree['loc'].iloc[r]]

		rDisp = random.uniform(0,1)
		disp_pool = []
		if rDisp < Pdisp:
			disp_pool = g.getNeighbours(old_pop.glob_index)
		else:
			disp_pool = g.global_grid[old_pop.glob_index].populations
			disp_pool = np.delete(disp_pool,old_pop.loc_index)

		if(TempSpeciation):
			Pspec_i = Pspec * tree['prob'].iloc[r] / boltz_min
		else:
			Pspec_i = Pspec

		rSpec = random.uniform(0,1)

		if(TempNiches):
			valid = False

			vcount = 0

			while(not valid):

				rNew = random.randint(0,len(disp_pool)-1)
				new_pop = disp_pool[rNew]

				T_old = 303.15 - (1/3)*np.absolute(g.global_grid[old_pop.glob_index].lat)
				T_new = 303.15 - (1/3)*np.absolute(g.global_grid[new_pop.glob_index].lat)

				temp_ancest_min = max(old_pop.T_min,new_pop.T_min) 
				temp_ancest_max = min(old_pop.T_max,new_pop.T_max)

				'''
				print('T_old ' + str(303.15 - (1/3)*np.absolute(g.global_grid[old_pop.glob_index].lat)))
				print('T_old_min ' + str(old_pop.T_min))
				print('T_old_max ' + str(old_pop.T_max))
				print('T_new_min ' + str(new_pop.T_min))
				print('T_new_max ' + str(new_pop.T_max))
				print('Temp_ancest_min ' + str(temp_ancest_min))
				print('Temp_ancest_max ' + str(temp_ancest_max))
				'''

				if(rSpec < Pspec_i):
					#print('Speciation')

					if(temp_ancest_max > temp_ancest_min and T_old < temp_ancest_max and T_old > temp_ancest_min and T_new < temp_ancest_max and T_new > temp_ancest_min):
						#print('valid')
						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_min = temp_ancest_min
						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_max = temp_ancest_max

						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_min = T_old - NicheWidth
						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_max = T_old + NicheWidth

						valid = True
						new_spec = Species(len(spec_list),old_pop.glob_index,old_pop.loc_index,T_old-NicheWidth,T_old+NicheWidth)
						spec_list.append(new_spec)
						IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index) & (IDlist['species']==-1), 'species'] = new_spec.order
						tree = tree.drop(tree.index[r])

				elif(len(tree[(tree['glob']==new_pop.glob_index) & (tree['loc']==new_pop.loc_index)]) != 0):
					#print('remove')
					'''
					print(temp_ancest_min)
					print(temp_ancest_max)
					print(T_old)
					print(T_old > temp_ancest_min)
					print(valid)
					'''
					if(temp_ancest_max > temp_ancest_min and T_old < temp_ancest_max and T_old > temp_ancest_min and T_new < temp_ancest_max and T_new > temp_ancest_min):

						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_min = T_old - NicheWidth
						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_max = T_old + NicheWidth

						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_min = temp_ancest_min
						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_max = temp_ancest_max

						glob_list = copy.deepcopy(IDlist['glob'])
						loc_list = copy.deepcopy(IDlist['loc'])
						IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'glob' ] = new_pop.glob_index
						IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'loc' ] = new_pop.loc_index
						tree = tree.drop(tree.index[r])

						valid = True

				else:
					#print('replace')
		
					if(temp_ancest_max > temp_ancest_min and T_old < temp_ancest_max and T_old > temp_ancest_min and T_new < temp_ancest_max and T_new > temp_ancest_min):

						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_min = T_old - NicheWidth
						g.global_grid[old_pop.glob_index].populations[old_pop.loc_index].T_max = T_old + NicheWidth

						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_min = temp_ancest_min
						g.global_grid[new_pop.glob_index].populations[new_pop.loc_index].T_max = temp_ancest_max

						valid = True
						glob_list = copy.deepcopy(IDlist['glob'])
						loc_list = copy.deepcopy(IDlist['loc'])
						IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'glob' ] = new_pop.glob_index
						IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'loc' ] = new_pop.loc_index
						tree.iloc[r,tree.columns.get_loc('glob')] = new_pop.glob_index
						tree.iloc[r,tree.columns.get_loc('loc')] = new_pop.loc_index
						if(TempTurnover):
							T = 303.15 - (1/3)*np.absolute(g.global_grid[new_pop.glob_index].lat)
							tree.iloc[r,tree.columns.get_loc('prob')] = np.exp(-0.65 / 8.617e-5 / T)
				
				#print(valid)
				if( (new_pop.T_min > T_old or new_pop.T_max < T_old) and valid == True):
					print('T_old ' + str(303.15 - (1/3)*np.absolute(g.global_grid[old_pop.glob_index].lat)))
					print('T_old_min ' + str(old_pop.T_min))
					print('T_old_max ' + str(old_pop.T_max))
					print('T_new_min ' + str(new_pop.T_min))
					print('T_new_max ' + str(new_pop.T_max))
					input()

				

				vcount += 1

				'''
				if(rSpec < Pspec_i or len(tree[(tree['glob']==new_pop.glob_index) & (tree['loc']==new_pop.loc_index)]) != 0):
					temp_min_new = IDlist['temp_min'].iloc[new_pop.glob_index*Nloc + new_pop.loc_index]
					temp_max_new = IDlist['temp_max'].iloc[new_pop.glob_index*Nloc + new_pop.loc_index]
					print('new ',str(temp_min_new) + '\t' + str(temp_max_new))

					temp_min_old = tree['temp_min'].iloc[r]
					temp_max_old = tree['temp_max'].iloc[r]
					print('old ',str(temp_min_old) + '\t' + str(temp_max_old))

					T = 303.15 - (1/3)*np.absolute(g.global_grid[tree['glob'].iloc[r]].lat) 
					print('T ',str(T))
					temp_min_new = min(temp_min_new,T)
					temp_max_new = max(temp_max_new,T)
					print('new ',str(temp_min_new) + '\t' + str(temp_max_new))

					if(rSpec < Pspec_i):
						print('speciation')
						temp_min_final = temp_min_new
						temp_max_final = temp_max_new
					else:
						print('remove')
						temp_min_final = min(temp_min_new,temp_min_old)
						temp_max_final = max(temp_max_new,temp_max_old)

					print('final ',str(temp_min_final), '\t', str(temp_max_final))

					temp_diff_final = temp_max_final - temp_min_final
					if(temp_diff_final > 2*NicheWidth):
						valid = False

					else:
						if(rSpec < Pspec_i):
							if(np.absolute(T - temp_min_final) <= NicheWidth and np.absolute(T - temp_max_final) <= NicheWidth):
								valid = True
								new_spec = Species(len(spec_list),old_pop.glob_index,old_pop.loc_index,T-NicheWidth,T+NicheWidth)
								spec_list.append(new_spec)
								IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index), 'temp_min'] = temp_min_new
								IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index), 'temp_max'] = temp_max_new
								IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index) & (IDlist['species']==-1), 'species'] = new_spec.order
								tree = tree.drop(tree.index[r])
							else:
								valid = False

						else:
							valid = True
							IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index), 'temp_min'] = temp_min_final
							IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index), 'temp_max'] = temp_max_final
							glob_list = copy.deepcopy(IDlist['glob'])
							loc_list = copy.deepcopy(IDlist['loc'])
							IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'glob' ] = new_pop.glob_index
							IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'loc' ] = new_pop.loc_index
							tree = tree.drop(tree.index[r])

				else:
					print('replace')
					temp_min_old = tree['temp_min'].iloc[r]
					temp_max_old = tree['temp_max'].iloc[r]

					print('old ',str(temp_min_old),'\t',str(temp_max_old))

					T = 303.15 - (1/3)*np.absolute(g.global_grid[new_pop.glob_index].lat)
					
					temp_min_final = min(temp_min_old,T)
					temp_max_final = max(temp_max_old,T)

					print('final ',str(temp_min_final),'\t',str(temp_max_final))

					if(temp_max_final - temp_min_final > 2*NicheWidth):
						valid = False
					else:
						valid = True
						IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index), 'temp_min'] = temp_min_final
						IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index), 'temp_max'] = temp_max_final
						glob_list = copy.deepcopy(IDlist['glob'])
						loc_list = copy.deepcopy(IDlist['loc'])
						IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'glob' ] = new_pop.glob_index
						IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'loc' ] = new_pop.loc_index
						tree.iloc[r,tree.columns.get_loc('glob')] = new_pop.glob_index
						tree.iloc[r,tree.columns.get_loc('loc')] = new_pop.loc_index
						tree.iloc[r,tree.columns.get_loc('temp_min')] = temp_min_final
						tree.iloc[r,tree.columns.get_loc('temp_max')] = temp_max_final
						if(TempTurnover):
							T = 303.15 - (1/3)*np.absolute(g.global_grid[new_pop.glob_index].lat)
							tree.iloc[r,tree.columns.get_loc('prob')] = np.exp(-0.65 / 8.617e-5 / T)

				vcount += 1
				if(vcount > 1000):
					for pop in disp_pool:
						print(str(IDlist['temp_min'].iloc[pop.glob_index*Nloc + pop.loc_index]),'\t',str(IDlist['temp_max'].iloc[pop.glob_index*Nloc + pop.loc_index]))
					input()
				'''
		else:
			if(rSpec < Pspec_i):
				print('speciation')
				print(len(spec_list))
				new_spec = Species(len(spec_list),old_pop.glob_index,old_pop.loc_index)
				spec_list.append(new_spec)
				IDlist.loc[ (IDlist['glob']==old_pop.glob_index) & (IDlist['loc']==old_pop.loc_index) & (IDlist['species']==-1), 'species'] = new_spec.order
				tree = tree.drop(tree.index[r])
			else:
				print('dispersal')
				rNew = random.randint(0,len(disp_pool)-1)
				new_pop = disp_pool[rNew]
				print('new')
				print(str(new_pop.glob_index) + '\t' + str(new_pop.loc_index))
				if len(tree[(tree['glob']==new_pop.glob_index) & (tree['loc']==new_pop.loc_index)]) != 0:
					print('remove')
					tree = tree.drop(tree.index[r])
				else:
					print('change')
					tree.iloc[r,tree.columns.get_loc('glob')] = new_pop.glob_index
					tree.iloc[r,tree.columns.get_loc('loc')] = new_pop.loc_index
					if(TempTurnover):
						T = 303.15 - (1/3)*np.absolute(g.global_grid[new_pop.glob_index].lat)
						tree.iloc[r,tree.columns.get_loc('prob')] = np.exp(-0.65 / 8.617e-5 / T)

				glob_list = copy.deepcopy(IDlist['glob'])
				loc_list = copy.deepcopy(IDlist['loc'])

				IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'glob' ] = new_pop.glob_index
				IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'loc' ] = new_pop.loc_index

		count += 1

	#print(tree)
	#print(IDlist)

	final_pop = g.global_grid[tree['glob'].iloc[0]].populations[tree['loc'].iloc[0]]
	final_spec = Species(len(spec_list),tree.iloc[0]['glob'],tree.iloc[0]['loc'])
	spec_list.append(final_spec)
	IDlist.loc[ (IDlist['glob']==tree['glob'].iloc[0]) & (IDlist['loc']==tree['loc'].iloc[0]) & (IDlist['species']==-1), 'species'] = final_spec.order

	#print(IDlist)


	c = 0
	for i in range(0,len(g.global_grid)):
		for j in range(0,len(g.global_grid[i].populations)):
			g.global_grid[i].populations[j].species = IDlist['species'][c]
			c += 1

	lon_list = np.zeros(len(g.global_grid))
	lat_list = np.zeros(len(g.global_grid))

	spec_array = np.zeros( (len(g.global_grid),len(spec_list)) )

	for i in range(0,len(lat_list)):
		lon_list[i] = g.global_grid[i].lon
		lat_list[i] = g.global_grid[i].lat

		for j in range(0,len(spec_list)):
			n = 0
			for x in range(0,len(g.global_grid[i].populations)):
				if(g.global_grid[i].populations[x].species == spec_list[j].order):
					n += 1	
			spec_array[i][j] = n

	ID_list = np.zeros(len(spec_list))
	for i in range(0,len(ID_list)):
		ID_list[i] = spec_list[i].order


	darray = np.concatenate((lon_list, lat_list)).reshape((-1, 2), order='F')
	darray = np.concatenate( (darray,spec_array) ,axis=1)
	names_out = ['lon','lat'] + ['spec_' + str(int(ID)) for ID in ID_list] 		
	df_out = pd.DataFrame(data=darray, columns=names_out)
			
	df_out.to_csv('Output/grid_' + str(int(k)).zfill(4) + '.csv')

t2 = time.time()

print(t2 - t1)