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
The code is executable, and detailed instructions on its use are available in the README files.  
The project contains multiple branches. Only a few are fully commented. Other branches are unfinished attempts not relevant for the thesis results.  

### Finished branches

#### master

The base code for the reconstructed Worm and Tittensor (2018) model.  
The code only contains the purely neutral death, birth, migration and speciation mechanics in a forward algorithm.  
The modules are fully commented, accompanied by a README file and are easily executable.  
No simulation run is added due to the size of the files, but a simulation can be run by the user.

#### dataRun

The code for a simulation of the empirical grids: NABB and CPR.  
The code shows the implementation of the additional physics, the coalescence algorithm and the parallelisation.  
The modules are fully commented, accompanied by a README file and are easily executable.  
An example of a simulation run is added in the 'OutputExample' directory.  
Additionally, a simulation can be run by the user.

#### theoryScaling

The code for a simulation of dummy grids for testing the Mechanistic Model.  
The code shows the implementation of the additional physics, the coalescence algorithm and the parallelisation.  
Additionally, the output modes are adapted, and extra plots are generated.  
The modules are fully commented, accompanied by a README file and are easily executable.  
An example of a simulation run is added in the 'OutputExample' directory.  
Additionally, a simulation can be run by the user.

### Unfinished branches

#### coalescence

First implementation of the coalescence algorithm.  
The code only contains the purely neutral death, birth, migration and speciation mechanics.  
The README file explains how to run a simulation although the code is not cleaned up.

#### coalescenceNewPhysics

Implementation of the additional physics in the coalescence algorithm.  
The code shows the implementation of: 

- Temperature-dependent death rates  
- Temperature-dependent speciation rates  
- Habitat area scaling of the local number of populations  
- Temperature ranges limiting migration 

The README file explains how to run a simulation although the code is not cleaned up.

#### newPhysics

Implementation of the additional physics in a forward algorithm.  
The code shows the implementation of: 

- Temperature-dependent death rates  
- Temperature-dependent speciation rates  
- Habitat area scaling of the local number of populations  
- Temperature ranges limiting migration  

The README file explains how to run a simulation although the code is not cleaned up.

#### parallelRun

Implementation of the parallelisation of the coalescence algorithm.  
The README file explains how to run a simulation although the code is not cleaned up.  

#### theoryScalingForward

Simulation of the dummy grid for testing the Mechanistic Model with a forward algorithm.  
The code is used for comparison with the output on the 'theoryScaling' branch.  
The README file explains how to run a simulation although the code is not cleaned up.

#### voter

Implementation of a Voter model in the same grid as the Worm and Tittensor (2018) model.  
The Voter Model is an alternative neutral model.  
The simulation is run with the coalescence algorithm.  
The README file explains how to run a simulation although the code is not cleaned up. 

---

## Instructions on the voter branch

Implementation of a Voter model in the same grid as the Worm and Tittensor (2018) model.  
The Voter Model is an alternative neutral model.  
The simulation is run with the coalescence algorithm.  
The code is not cleaned up but can be run as follows:  

1. Download the files and put them in a local directory.
2. Execute the 'tree.py' python file: $ python3 tree.py .  
!!! The code uses all available processor cores. Make sure that your computer can run this.
3. Change parameters in 'parameters.py' as you please.
4. Output is generated in the 'Output' directory.
5. Plots are generated in the 'plot.R' script. This has to be done manually.