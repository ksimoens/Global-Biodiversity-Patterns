from parameters import*
from Grid import*

import numpy as np
import pandas as pd
import random
import copy

g = Grid()
g.fillGrid()

class Species():

	def __init__(self,o,i_glo,i_loc):
		self.order = o
		self.glob_index = i_glo
		self.loc_index = i_loc

spec_list = []

glob_index = []
loc_index = []

for i in range(0,len(g.global_grid)):
	for j in range(0,len(g.global_grid[i].populations)):
		glob_index.append(i)
		loc_index.append(j)

tree = pd.DataFrame(data=np.vstack( (glob_index,loc_index)).T,columns=['glob','loc'])
print(tree)
IDlist = pd.DataFrame(data=np.vstack( (glob_index,loc_index)).T,columns=['glob','loc'])
IDlist['species'] = -1
print(IDlist)

count = 0

while(len(tree) > 1):

	print("iteration:" + '\t' + str(count))
	r = random.randint(0,len(tree)-1)
	old_pop = g.global_grid[tree.iloc[r]['glob']].populations[tree.iloc[r]['loc']]
	print('old')
	print(tree.iloc[r])

	rDisp = random.uniform(0,1)
	disp_pool = []
	if rDisp < Pdisp:
		disp_pool = g.getNeighbours(old_pop.glob_index)
	else:
		disp_pool = g.global_grid[old_pop.glob_index].populations
		disp_pool = np.delete(disp_pool,old_pop.loc_index)

	rSpec = random.uniform(0,1)
	if(rSpec < Pspec):
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
			tree.iloc[r]['glob'] = new_pop.glob_index
			tree.iloc[r]['loc'] = new_pop.loc_index

		glob_list = copy.deepcopy(IDlist['glob'])
		loc_list = copy.deepcopy(IDlist['loc'])

		IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'glob' ] = new_pop.glob_index
		IDlist.loc[ (glob_list==old_pop.glob_index) & (loc_list==old_pop.loc_index),'loc' ] = new_pop.loc_index


	count += 1

print(tree)
print(IDlist)

final_pop = g.global_grid[tree.iloc[0]['glob']].populations[tree.iloc[0]['loc']]
final_spec = Species(len(spec_list),tree.iloc[0]['glob'],tree.iloc[0]['loc'])
spec_list.append(final_spec)
IDlist.loc[ (IDlist['glob']==tree.iloc[0]['glob']) & (IDlist['loc']==tree.iloc[0]['loc']) & (IDlist['species']==-1), 'species'] = final_spec.order

print(IDlist)


c = 0
for i in range(0,len(g.global_grid)):
	for j in range(0,len(g.global_grid[i].populations)):
		g.global_grid[i].populations[j].species = IDlist['species'][c]
		c += 1

for i in range(0,len(g.global_grid)):
	for j in range(0,len(g.global_grid[i].populations)):
		print(str(i) + "\t" + str(j) + "\t" +  str(g.global_grid[i].populations[j].species))
		

