import numpy as np

Nlon = 45
Nlat = 14

Nloc = 16
Nspec0 = 1

Nlon = int(Nlon*np.sqrt(Nloc))
Nlat = int(Nlat*np.sqrt(Nloc))

Pdisp = 0.1
Pspec = 0.0001