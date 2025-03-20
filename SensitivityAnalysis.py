"""
Author: Joe Zavorskas
Date Started: 03/10/2025
Last Edit: 03/19/2025

This file will perform sensitivity analysis on the two genome-scale metabolic models (GSMMs) used in
this project. The diatom-only and symbiotic diatom/cyanobacteria models will be run through
flux variabilitity analysis (FVA) under multiple nutrient conditions, where relevant:

    1) Rubisco (Carbon) limitation
    2) Nitrogen limitation
    3) Silicon limitation
    4) Iron limitiation

Iron limitation is not relevant to the diatom-only model, as the diatom GSMM does not consider the
iron uptake required to sustain photosynthesis.

Each of these flux variability analyses will be run in day and night conditions.

"""

# Import necessary packages:

import cobra
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
from datetime import datetime, timedelta
from scipy.interpolate import RegularGridInterpolator
import netCDF4
from IPython.display import clear_output
import sys
from time import time

sys.path.append("C:\\Users\\josep\\OneDrive - University of Connecticut\\GRADUATE SCHOOL\\Srivastava Group\\~Organized Grad School\\Ch.6 dFBA Model\\VSCode Port\\Code_NoDepth")
sys.path.append("C:\\Users\\josep\\OneDrive - University of Connecticut\\GRADUATE SCHOOL\\Srivastava Group\\~Organized Grad School\\Ch.6 dFBA Model\\VSCode Port\\Code_NoDepth\\Auxiliary_Functions")

from Auxiliary_Functions import Light, Extras, Input, LoadModel, MergeLoadModel

# Initialize the inputs dictionary which will be used through dFBA simulation.
[Inputs, start_date] = Input.Input()

[ML, DML, refmodel, IYSi, INSi, Inputs] = MergeLoadModel.MergeLoadModel(Inputs)

Scenario = ["RuBisCo","Nitrogen","Silicon","Iron"]
LightType = ["Day","Night"]

NitrogenDay = np.array([[9,0.25,9,9],[9,0.25,9,9]])
SiliconDay = np.array([[4,4,0.25,4],[4,4,0.25,4]])
IronDay = np.array([[0.004,0.004,0.004,0.0001],[0.004,0.004,0.004,0.0001]])
T = 5 # degrees Celsius

# Calculate the saturation concentration of carbon dioxide (in mM).
Henry = Extras.henry(T)
CO2sat = Henry[0]*(4.21e-4) # Partial pressure of CO2 in atm. # convert to uM
O2sat = Henry[1]*(.2105) # convert to uM
N2sat = Henry[2]*(.781)
n = 1.9

starttime = time()
Iphoto = Inputs['Iphoto']

for idy, Nutrient in enumerate(Scenario):

    for idx, Type in enumerate(LightType):

        print("Nutrient:", Nutrient)
        print("Type:", Type)

        if idx == 0:

            for k in Iphoto.keys():
                DML.reactions.get_by_id(k).upper_bound = 0
                DML.reactions.get_by_id(k).lower_bound = -1
                DML.reactions.get_by_id(k).upper_bound = -0.9999
                ML.reactions.get_by_id(k).upper_bound = 0
                ML.reactions.get_by_id(k).lower_bound = -1
                ML.reactions.get_by_id(k).upper_bound = -0.9999

            ML.reactions.EX_co2_e.upper_bound = 0
            ML.reactions.EX_co2_e.lower_bound = -Extras.RLPD(CO2sat*4,O2sat,T)
            ML.reactions.EX_chryso_e.lower_bound = 0
            DML.reactions.EX_chryso_e.lower_bound = 0
            DML.reactions.EX_co2_e.upper_bound = 0
            DML.reactions.EX_co2_e.lower_bound = -Extras.RLPD(CO2sat*4,O2sat,T)
        
        else:

            for k in Iphoto.keys():
                DML.reactions.get_by_id(k).upper_bound = 0
                DML.reactions.get_by_id(k).lower_bound = Iphoto[k]*-1.
                DML.reactions.get_by_id(k).upper_bound = Iphoto[k]*-0.9999
                ML.reactions.get_by_id(k).upper_bound = 0
                ML.reactions.get_by_id(k).lower_bound = Iphoto[k]*-1.
                ML.reactions.get_by_id(k).upper_bound = Iphoto[k]*-0.9999

            DML.reactions.EX_chryso_e.lower_bound = -0.5
            ML.reactions.EX_chryso_e.lower_bound = -0.5
            DML.reactions.EX_co2_e.lower_bound = 0
            ML.reactions.EX_co2_e.lower_bound = 0
            ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound = 0

        ML.reactions.EX_no3_e.lower_bound = -((1.474*NitrogenDay[idx,idy])/(6.14+NitrogenDay[idx,idy]))
        DML.reactions.EX_no3_e.lower_bound = -((1.474*NitrogenDay[idx,idy])/(6.14+NitrogenDay[idx,idy]))

        ML.reactions.EX_sio4h4_e.lower_bound = -((1.961*(SiliconDay[idx,idy]**n))/((8.1**n)+(SiliconDay[idx,idy]**n)))
        DML.reactions.EX_sio4h4_e.lower_bound = -((1.961*(SiliconDay[idx,idy]**n))/((8.1**n)+(SiliconDay[idx,idy]**n)))

        ML.reactions.get_by_id("EX_cpd10515[e]").lower_bound = -((0.0024*IronDay[idx,idy])/(0.0023+IronDay[idx,idy]))

        DML.optimize()
        ML.optimize()

        print("Diatom-Only")
        FVADML = cobra.flux_analysis.variability.flux_variability_analysis(DML)
        print("Symbiotic")
        FVAML = cobra.flux_analysis.variability.flux_variability_analysis(ML)

        DMLFileName = "./SensitivityOutputs/DiatomOnly_FVA_"+Nutrient+"Limit_"+Type+".csv"
        MLFileName = "./SensitivityOutputs/Symbiotic_FVA_"+Nutrient+"Limit_"+Type+".csv"

        FVADML.to_csv(DMLFileName, index=True) 
        FVAML.to_csv(MLFileName, index=True) 

        print("Completion:", (((idy*4+idx)+1)/8)*100,"%")
        print("Elapsed Time:", (time() - starttime))

        