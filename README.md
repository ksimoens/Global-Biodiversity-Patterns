# Global Biodiversity Patterns
## Master Thesis IMBRSea

The Physics of Biodiversity:
exploring the dynamics behind spatial biodiversity patterns
<br/>  
MECHANISTIC MODEL SIMULATIONS
<br/>  
Student: Kobe Simoens  
Contact: kobe.simoens@imbrsea.eu

### Promotors

Sandro Azaele  
Università di Padova  
https://orcid.org/0000-0002-5153-4833  
Contact: sandro.azaele@unipd.it  
<br/>
Jan Ryckebusch  
Universiteit Gent  
https://orcid.org/0000-0001-7750-1522  
Contact: jan.ryckebusch@ugent.be

### Supervisors

Samir Suweis  
Università di Padova  
https://orcid.org/0000-0002-1603-8375  
Contact: samir.suweis@unipd.it  
<br/>
Jan Ryckebusch  
Universiteit Gent  
https://orcid.org/0000-0001-7750-1522  
Contact: jan.ryckebusch@ugent.be  
<br/>
Éric Goberville  
Sorbonne Université  
https://orcid.org/0000-0002-1843-7855  
Contact: eric.goberville@upmc.fr  
<br/>
Derek Tittensor  
Dalhousie University   
https://orcid.org/0000-0002-9550-3123  
Contact: derek.tittensor@dal.ca  

---

## Instructions

This GitHub repository contains all information on the simulations of the Mechanistic Model of Worm and Tittensor (2018).  
The code is executable and detailed instructions on its use is available in the README files.  
The project contains multiple branches. Only a few are fully commented. Other branches are unfinished attempts not relevant for the thesis results.  

### Finished branches

#### Master

The base code for the reconstructed Worm and Tittensor (2018) model.  
The code only contains the purely neutral death, birth, migration and speciation mechanics.  
The modules are fully commented, accompanied by a README file and easily executable.  
No simulation run is added due to the size of the files, but a simulation can be run by the user.

### Unfinished branches

---

## Instructions on the Master branch

This branch contains the code for the base model.  
The code is easily executable:  

1. Download the files and put them in a local directory  
2. Execute the 'runSimulation.sh' bash file: $ ./runSimulation.sh  
The full run will take around 40 minutes.  
3. Change parameters in 'parameters.py' as you please.  
However, there is no guarantee that all combinations of parameters will run equally efficiently.  
4. Output is generated in the 'Output' directory; including some plots in the 'plots' subdirectory.

### Files

- **\_\_pycache\_\_**: directory for the pycache files, which ensure communication between the different modules.  
These files are not strictly necessary and will be generated automatically when running the code.

- **Grid.py**: module that contains the 'Grid' class.  
An object of the Grid class instantiates the simulation grid.  
The class contains information on the grid geometry and size.  
The grid also contains communities of populations.  
The Grid class holds the method for a turnover, which is the working horse of the simulation.  
A Grid object is therefore updated each turnover and can be printed out as simulation output.

- **Local.py**: module that contains the 'Local' class.  
An object of the Local class instantiates a local community of populations / a grid cell.  
The Local object contains populations represented by integers, which indicate the species identity.

- **land_all_data_inc_working.csv**: definition of the simulation grid by Worm and Tittensor (2018).  
The file is mainly used for reading in the grid geometry including the cell coordinates.  
For details on the data format, the user is referred to the authors of the data file.

- **main.py**: runs the actual simulation.  
Execute this file to run a simulation.  

- **parameters.py**: contains all relevant simulation parameters.  
Grid parameters, turnover parameters and numerical parameters.  
The parameters are written to a txt file after each simulation for future reference.

- **plot.R**: R script to make some plots based on the model output.  
The script automatically reads the model output.  
Some plot parameters might not be suited for other simulation parameter values.

- **runSimulation.sh**: simple executable which runs the simulation and makes the plots.  
Execute this file or run main.py and plot.R manually.

- **Output**: output format:  

	- **grid_XXXX.csv**: csv file with a diversity matrix for the simulation grid.  
	XXXX represents a unique identifier and can be the simulation time at which the grid was printed.  
	Columns are:
		- *lon*: longitudinal coordinate of the grid cell as defined by Worm and Tittensor (2018)  
		- *lat*: latitudinal coordinate of the grid cell as defined by Worm and Tittensor (2018)  
		- *$spec$*:  
		All remaining columns represent a unique species.  
		Column names are the unique simulation identifiers for the species.  
		Values are the total number of entries of a particular species found in a particular grid cell.  

	- **PARAM_file.txt**: a txt file with information on the simulation parameters.
