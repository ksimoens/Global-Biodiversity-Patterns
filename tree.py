from parameters import*
from Grid import*

import numpy as np
import pandas as pd
import random
import copy
import time

g = Grid()
g.fillGrid()



class Species():

	def __init__(self,o,i_glo,i_loc):
		self.order = o
		self.glob_index = i_glo
		self.loc_index = i_loc

spec_list = []

lat_index = []
lon_index = []


for i in range(0,len(g.global_grid)):
	for j in range(0,len(g.global_grid[i])):
		lat_index.append(i)
		lon_index.append(j)



tree = pd.DataFrame(data=np.vstack( (lat_index,lon_index)).T,columns=['lat','lon'])
print(tree)
IDlist = pd.DataFrame(data=np.vstack( (lat_index,lon_index)).T,columns=['lat','lon'])
IDlist['species'] = -1
print(IDlist)

count = 0

t1 = time.time()

while(len(tree) > 1):

	print("iteration:" + '\t' + str(count) + '\t' + str(len(tree)))
	r = random.randint(0,len(tree)-1)
	i_lat = tree.iloc[r]['lat']
	i_lon = tree.iloc[r]['lon']
	old_pop = g.global_grid[i_lat][i_lon]
	print('old')
	print(tree.iloc[r])

	rDisp = random.uniform(0,1)
	disp_pool = g.getNeighbours(i_lat,i_lon)
	'''
	if rDisp < Pdisp:
		disp_pool = g.getNeighbours(old_pop.glob_index)
	else:
		disp_pool = g.global_grid[old_pop.glob_index].populations
		disp_pool = np.delete(disp_pool,old_pop.loc_index)
	'''
	rNew = random.randint(0,len(disp_pool)-1)
	new_pop = disp_pool[rNew]
	print('new')
	print(str(new_pop.lat_index) + '\t' + str(new_pop.lon_index))


	rSpec = random.uniform(0,1)
	if(rSpec < Pspec):
		print('speciation')
		print(len(spec_list))
		new_spec = Species(len(spec_list),old_pop.lat_index,old_pop.lon_index)
		spec_list.append(new_spec)
		IDlist.loc[ (IDlist['lat']==old_pop.lat_index) & (IDlist['lon']==old_pop.lon_index) & (IDlist['species']==-1), 'species'] = new_spec.order
		tree = tree.drop(tree.index[r])
	else:
		print('dispersal')
		if len(tree[(tree['lat']==new_pop.lat_index) & (tree['lon']==new_pop.lon_index)]) != 0:
			print('remove')
			tree = tree.drop(tree.index[r])
		else:
			print('change')
			tree.iloc[r]['lat'] = new_pop.lat_index
			tree.iloc[r]['lon'] = new_pop.lon_index

	glob_list = copy.deepcopy(IDlist['lat'])
	loc_list = copy.deepcopy(IDlist['lon'])

	IDlist.loc[ (glob_list==old_pop.lat_index) & (loc_list==old_pop.lon_index),'lat' ] = new_pop.lat_index
	IDlist.loc[ (glob_list==old_pop.lat_index) & (loc_list==old_pop.lon_index),'lon' ] = new_pop.lon_index


	count += 1

print(tree)
print(IDlist)

final_pop = g.global_grid[tree.iloc[0]['lat']][tree.iloc[0]['lon']]
final_spec = Species(len(spec_list),tree.iloc[0]['lat'],tree.iloc[0]['lon'])
spec_list.append(final_spec)
IDlist.loc[ (IDlist['lat']==tree.iloc[0]['lat']) & (IDlist['lon']==tree.iloc[0]['lon']) & (IDlist['species']==-1), 'species'] = final_spec.order

print(IDlist)


'''
dummy_spec = np.zeros(Nlon*Nlat)
for i in range(0,len(dummy_spec)):
	dummy_spec[i] = random.randint(0,500)

uniq_list = np.unique(dummy_spec,return_counts=True)

spec_list = []

for i in range(0,len(uniq_list)):
	spec_list.append(Species(uniq_list[i],0,0))
'''

c = 0
for i in range(0,len(g.global_grid)):
	for j in range(0,len(g.global_grid[i])):
		g.global_grid[i][j].species = IDlist['species'][c]
		c += 1

nlon = int(Nlon/np.sqrt(Nloc))
nlat = int(Nlat/np.sqrt(Nloc))

df_grid = pd.read_csv('land_all_data_inc_working.csv').iloc[:,[0,5,6]]
mapped_grid = np.empty([nlat,nlon]).tolist()

c = 0
for i in range(0,len(mapped_grid)):
	for j in range(0,len(mapped_grid[i])):
		mapped_grid[i][j] = Local(df_grid.iloc[c,1],df_grid.iloc[c,2])

		for k in range(i*int(np.sqrt(Nloc)),i*int(np.sqrt(Nloc))+int(np.sqrt(Nloc))):
			for l in range(j*int(np.sqrt(Nloc)),j*int(np.sqrt(Nloc))+int(np.sqrt(Nloc))):
				mapped_grid[i][j].populations.append(g.global_grid[k][l])
		c += 1

lon_list = np.zeros(len(mapped_grid)*len(mapped_grid[0]))
lat_list = np.zeros(len(mapped_grid)*len(mapped_grid[0]))

spec_array = np.zeros( (len(lon_list),len(spec_list)) )

c = 0
for i in range(0,len(mapped_grid)):
	for j in range(0,len(mapped_grid[i])):
		lon_list[c] = mapped_grid[i][j].lon
		lat_list[c] = mapped_grid[i][j].lat

		for k in range(0,len(spec_list)):
			n = 0
			for x in range(0,len(mapped_grid[i][j].populations)):
				if(mapped_grid[i][j].populations[x].species == spec_list[k].order):
					n += 1	
			

			spec_array[c][k] = n

		c += 1

ID_list = np.zeros(len(spec_list))
for i in range(0,len(ID_list)):
	ID_list[i] = spec_list[i].order


darray = np.concatenate((lon_list, lat_list)).reshape((-1, 2), order='F')
darray = np.concatenate( (darray,spec_array) ,axis=1)
names_out = ['lon','lat'] + ['spec_' + str(int(ID)) for ID in ID_list] 		
df_out = pd.DataFrame(data=darray, columns=names_out)
		
df_out.to_csv('grid_' + str(int(0)).zfill(4) + '.csv')

t2 = time.time()

print(t2 - t1)
		

