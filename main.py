# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Simulation of Mechanistic Model
#####################################
# MAIN EXECUTION of the simulation
#####################################


# ------------ IMPORT FROM OTHER FILES ---------------------

from tree import*
from parameters import*

# ----------------------------------------------------------


# ----------------- IMPORT MODULES -------------------------

import time
import os
import shutil
import multiprocessing as mp

# ----------------------------------------------------------



# create a directory for output files
if(os.path.exists("Output")):
	shutil.rmtree("Output")
os.mkdir("Output")

# time the complete run
t1 = time.time()

# create a pool of child processes
# NCPU defined in parameters.py
pool = mp.Pool(NCPU)
# divide the replicate runs over all child cores
# replicateRun defined in tree.py
# Nrep defined in parameters.py
pool.map(replicateRun, [k for k in range(0,Nrep)])
# close the pool
pool.close()

t2 = time.time()
print('\nwall time: ' + str(round((t2-t1)/60,4)) + ' minutes')
print_parameters(t2-t1)
