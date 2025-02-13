"""
Author: Joe Zavorskas
Date Started: 6/26/2023
Last Edit: 04/11/2024

Main file for diatom/cyanobacteria simulations. Simulates the symbiosis between Chaetoceros and 
a cyanobacterium, which relieves nitrogen stress on Chaetoceros after nitrate has been depleted.

A detailed input instance is passed throughout the differential equation function defined below, which 
contains information required to calculate light intensity, and saved LP solutions for the diatom-only FBA
and the symbiotic FBA problems.

Please see "Symbiosis_dydt" for a more thorough understanding of this file.
"""
"""
y: 
0  [   C.socialis Biomass  ] mg dry weight/L
1  [    Cyano Biomass      ] mg dry weight/L
2  [   C.socialis Chryso   ] umol/L
3  [       Nitrate         ] umol/L
4  [       Silicon         ] umol/L
5  [        Iron           ] umol/L
6  [    Total Biomass      ] mg dry weight/L
7  [      Total CO2        ] umol/L
8  [      Total O2         ] umol/L
9  [      Total N2         ] umol/L
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

sys.path.append("C:\\Users\\josep\\OneDrive - University of Connecticut\\GRADUATE SCHOOL\\Srivastava Group\\~Organized Grad School\\Ch.6 dFBA Model\\VSCode Port\\Code_NoDepth")
sys.path.append("C:\\Users\\josep\\OneDrive - University of Connecticut\\GRADUATE SCHOOL\\Srivastava Group\\~Organized Grad School\\Ch.6 dFBA Model\\VSCode Port\\Code_NoDepth\\Auxiliary_Functions")
sys.path.append("C:\\Users\\josep\\OneDrive - University of Connecticut\\GRADUATE SCHOOL\\Srivastava Group\\~Organized Grad School\\Ch.6 dFBA Model\\VSCode Port\\Code_NoDepth\\Hierarchical_Optimization")
sys.path.append("C:\\Users\\josep\\OneDrive - University of Connecticut\\GRADUATE SCHOOL\\Srivastava Group\\~Organized Grad School\\Ch.6 dFBA Model\\VSCode Port\\Code_NoDepth\\Validation")

# Author-written functions:
from Code_NoDepth.Hierarchical_Optimization import MergefluxCalc, MergeDaySolve, MergeNightSolve
from Code_NoDepth.Hierarchical_Optimization import fluxCalc, DaySolve, NightSolve
from Code_NoDepth.Auxiliary_Functions import Light, Extras, Input, LoadModel, MergeLoadModel
from Symbiosis_dydt import Symbiosis_dydt as dydt

# Initialize the inputs dictionary which will be used through dFBA simulation.
[Inputs, start_date] = Input.Input()

#__________________________
### Main Simulation Section
#__________________________


lats = [55, 65, 75, 85]
#lats = [70]
time = []
latmap = []
Dbiomass = []
Cbiomass = []
NO3 = []
SiO4H4 = []
totalBM = []
totalCO2 = []
totalO2 = []
totalN2 = []

for idx,lat in enumerate(lats):

    print("Latitude: " + str(lat))
    global ML
    # Import three models, and set their constraints to default:
    # 1) diatom only, 2) symbiotic, 3) cyanobacteria only, for reference.
    [ML, DML, refmodel, IYSi, INSi, Inputs] = MergeLoadModel.MergeLoadModel(Inputs)

    # Store the interpolators inside of the Inputs dictionary.
    Inputs['IYSi'] = IYSi
    Inputs['INSi'] = INSi

    # Set initial values and flags. These assume that the simulation will start at daytime.
    Inputs['daycounter'] = 0
    Inputs['StartBool'] = True
    Inputs['FirstDay'] = True
    Inputs['FirstNight'] = False
    Inputs['FirstLight'] = False
    # Extra flag to trigger change between day and night at low light intensity.
    Inputs['LightLimit'] = 5
    # Extra flag to force 500% of basal nitrogen fixation to be transferred from cyanobacteria
    # to the diatom.
    Inputs['NH4Constraint'] = 5

    # Load sea-surface temperature values between March 1st and August 1st.
    # This array is at 2.5 degree latitude resolution, so ID is multiplied by 4.    
    temp = np.loadtxt('.\\Data\\sstdailyMarAug2020.csv',delimiter=',')[:,idx*4]
    Henry = Extras.henry(int(temp[0]))
    CO2sat = Henry[0]*(4.21e-4)
    O2sat = Henry[1]*(.2105)
    Inputs['Latitude'] = lat

    # Set the initial constraints of the reference cyanobacteria-only model. This model is
    # important because many of the transfer reactions are scaled off of its basal metabolism.
    refmodel.reactions.get_by_id("EX_cpd00011[e]").lower_bound = -Extras.RLPC(CO2sat*4,O2sat,5)
    refmodel.reactions.get_by_id("EX_cpd00242[e]").lower_bound = 0
    refmodel.reactions.get_by_id("EX_cpd00027[e]").lower_bound = 0
    refmodel.reactions.biomass_eq_33047__hc.lower_bound = 0
    
    print("Entering ODE Calculator")
    IVPSol = scipy.integrate.solve_ivp(dydt,[0, 4300],[0.0005, 0.0005, 0.055, 27, 5, 9, 0.001, 0, 0, 0],max_step=1,args=[Inputs,start_date,temp,ML,refmodel])

    time.append(IVPSol['t'])
    
    latinput = np.zeros((1,len(IVPSol['t'])))
    latinput.fill(lat)
    latmap.append(latinput)
    
    Dbiomass.append(IVPSol['y'][0])
    Cbiomass.append(IVPSol['y'][1])
    NO3.append(IVPSol['y'][4])
    SiO4H4.append(IVPSol['y'][5])
    totalBM.append(IVPSol['y'][6])
    totalCO2.append(IVPSol['y'][7])
    totalO2.append(IVPSol['y'][8])
    totalN2.append(IVPSol['y'][9])
    
timeflat = Extras.listarraytoarray(time)
latflat = Extras.flattenlat(latmap)
Dbioflat = Extras.listarraytoarray(Dbiomass)
Cbioflat = Extras.listarraytoarray(Cbiomass)
NO3flat = Extras.listarraytoarray(NO3)
SiO4H4flat = Extras.listarraytoarray(SiO4H4)
totalBMflat = Extras.listarraytoarray(totalBM)
totalCO2flat = Extras.listarraytoarray(totalCO2)
totalO2flat = Extras.listarraytoarray(totalO2)
totalN2flat = Extras.listarraytoarray(totalN2)

output = np.array([timeflat, latflat, Dbioflat, Cbioflat, NO3flat, SiO4H4flat, totalBMflat, totalCO2flat, totalO2flat, totalN2flat])
np.savetxt("SymbiosisOnly.csv",output,delimiter=',')

### Section Break Here

March = np.loadtxt("SymbiosisOnly.csv",delimiter=",")
lengths = []
#lats = [65, 67.5, 70, 72.5, 75, 77.5, 80, 82.5, 85]
lats = [55, 65, 75, 85]
#lats = [75, 85]
#lats = [70]

numrows = len(lats)

for lat in lats:
    lengths.append(len(np.where(March[1] == lat)[0]))

arraydict = {}
arraydict['arraytime'] = np.zeros((numrows,max(lengths)))
arraydict['arraylat'] = np.zeros((numrows,max(lengths)))
arraydict['arrayDbio'] = np.zeros((numrows,max(lengths)))
arraydict['arrayCbio'] = np.zeros((numrows,max(lengths)))
arraydict['arrayNO3'] = np.zeros((numrows,max(lengths)))
arraydict['arraySiO4H4'] = np.zeros((numrows,max(lengths)))
arraydict['arraytotalBM'] =  np.zeros((numrows,max(lengths)))
arraydict['arraytotalCO2'] =  np.zeros((numrows,max(lengths)))
arraydict['arraytotalO2'] =  np.zeros((numrows,max(lengths)))
arraydict['arraytotalN2'] =  np.zeros((numrows,max(lengths)))
count = 0

for entry in arraydict.keys():
    
    for idx, lat in enumerate(lats):
        ids = np.where(March[1] == lat)
        hold = March[count][ids[0][0]:ids[0][-1]]
        holdpad = np.pad(hold,(0,max(lengths)-len(ids[0])+1),'constant',constant_values=(0,hold[-1]))
        arraydict[entry][idx,:] = holdpad
        
    count += 1

### Section Break Here

x = arraydict['arraytime']
lats = [55,65,75,85]
variables = ['arrayDbio','arrayCbio','arrayNO3','arraySiO4H4','arraytotalBM','arraytotalCO2','arraytotalO2','arraytotalN2']
labels = ["Active Diatom Biomass \n (mgDW/L)", "Active Cyanobacteria Biomass \n (mgDW/L)", "Nitrate Concentration \n (umol/L)",
           "Silicic Acid Concentration \n (umol/L)","Total Biomass Generated \n (mgDW/L)", "Total CO2 Consumed \n (umol/L)",
           "Total O2 Produced \n (umol/L)", "Total N2 Consumed \n (umol/L)"]
filenames = ["DiatomBio","CyanoBio","Nitrate","Silicic","TotalBio","TotalCO2","TotalO2","TotalN2"]


for idx, name in enumerate(filenames):
    plt.figure()
            
    for idy, lat in enumerate(lats):
        plt.scatter(x[idy,:],arraydict[variables[idx]][idy,:],s=0.75)

    plt.xlabel("Time (hours)")
    plt.ylabel(labels[idx]) 

    plt.rc('font', size=16)          # controls default text sizes
    plt.rc('axes', titlesize=18)     # fontsize of the axes title
    plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
    plt.rc('legend', fontsize=16)    # legend fontsize
    plt.rc('figure', titlesize=12)  # fontsize of the figure title
    plt.savefig(("March-AugustNewBigNoLegend " + name + ".pdf"),dpi=1200,format='pdf')

### Section Break Here



