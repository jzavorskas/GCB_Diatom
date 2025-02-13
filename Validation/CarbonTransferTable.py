"""
Author: Joe Zavorskas
Date Started: 04/11/2024
Last Edit: 02/11/2025
_______________________________________________________________________
Script: Carbon Transfer Lookup Table (CarbonTransferTable) generation

This file, used in tandem with "CarbonTransferTable.py," generates a lookup table of the 
carbon transfer value that leads to equal growth rates for diatom and cyanobacteria at a 
range of temperatures, nitrate concentrations, and silicic acid concentrations.

See "CarbonTransferAbs.py" for further documentation.
"""

import CarbonTransferAbs, Extras, MergeLoadModel
import numpy as np
from IPython.display import clear_output
import Input

# Initialize inputs array and perform the first light intensity calculation.
[Inputs,start_date] = Input()
# Load the symbiotic model, adding transfer and storage reactions.
[ML, DML, refmodel, IYSi, INSi, Inputs] = MergeLoadModel(Inputs)

# Generate grids for temperature and concentrations over which to generate the lookup table.
# Nitrate and silicon concentration ranges were selected based on the concentrations each
# nutrient is known to be limiting (CO2 consumption is less than expected in rubisco-limited conditions)
Tn = np.linspace(-3,26,30)
NO3 = np.linspace(0,.5,3)
Sili = np.linspace(0.2,2,18)

# Create a 3-D array that will hold the calculated carbon transfer values as each variable changes.
datamin = np.zeros((len(Tn),len(NO3),len(Sili)))

for idx, temp in enumerate(Tn):     
    for idy, NO3ele in enumerate(NO3):
        for idz, Siliele in enumerate(Sili):
            
            # Given the current loop values for each variable, calculate the expected carbon transfer.
            sol = CarbonTransferAbs(temp,NO3ele,Siliele,ML,refmodel,Inputs)
    
            print("Temperature:", Tn[idx])
            print("NO3:", idy)
            print("Si(OH)4:", idz)
            print(sol)
            clear_output(wait=True)
    
            datamin[idx,idy,idz] = sol

# Reshape the 3-D array into 2-D so it can be saved as a .csv file.     
datamin_re = datamin.reshape(datamin.shape[0], -1)
# Save it
np.savetxt("CarbonTransfer.csv",datamin_re,delimiter=',')