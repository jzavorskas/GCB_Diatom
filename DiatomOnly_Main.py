"""
Author: Joe Zavorskas
Date Started: 6/26/2023
Last Edit: 04/11/2024

This file is the main file for a dynamic flux balance analysis (dFBA) simulation of diatom growth and nutrient consumption
in the Arctic. A flux balanced analysis linear program reports consumption or generation of metabolites, which are used to update nutrient concentrations in a system of ordinary differential equations (ODEs).

We have written a number of auxiliary files containing solutions to common dFBA problems such as:
    1) all fluxes other than objective function are non-unique, which can cause instability in ODEs. Solved via hierarchical optimization.
        (fluxCalc)
    2) dFBA often performs redundant optimizations. Solved by saving and scaling solutions. 
        (daySolve, nightSolve).

"""

# Import necessary packages:
import cobra
from cobra.io import save_matlab_model,write_sbml_model
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
from datetime import datetime, timedelta
from scipy.interpolate import RegularGridInterpolator
import netCDF4
from IPython.display import clear_output
from DiatomOnly_dydt import DiatomOnly_dydt as dydt
import sys

# Author-written functions:
from Code_NoDepth.Hierarchical_Optimization import fluxCalc, DaySolve, NightSolve
from Code_NoDepth.Auxiliary_Functions import Light, Extras, Input, LoadModel, Input

sys.path.append("C:\\Users\\josep\\OneDrive - University of Connecticut\\GRADUATE SCHOOL\\Srivastava Group\\~Organized Grad School\\Ch.6 dFBA Model\\VSCode Port\\Code_NoDepth")
sys.path.append("C:\\Users\\josep\\OneDrive - University of Connecticut\\GRADUATE SCHOOL\\Srivastava Group\\~Organized Grad School\\Ch.6 dFBA Model\\VSCode Port\\Code_NoDepth\\Auxiliary_Functions")
sys.path.append("C:\\Users\\josep\\OneDrive - University of Connecticut\\GRADUATE SCHOOL\\Srivastava Group\\~Organized Grad School\\Ch.6 dFBA Model\\VSCode Port\\Code_NoDepth\\Hierarchical_Optimization")
sys.path.append("C:\\Users\\josep\\OneDrive - University of Connecticut\\GRADUATE SCHOOL\\Srivastava Group\\~Organized Grad School\\Ch.6 dFBA Model\\VSCode Port\\Code_NoDepth\\Validation")

# Initialize the inputs dictionary which will be used through dFBA simulation.
[Inputs, start_date] = Input.Input()

#__________________________
### Main Simulation Section
#__________________________

global daycounter

# Initialize all arrays that will hold dFBA data.
lats = [55,65,75,85]
time = []
latmap = []
biomass = []
totalBM = []
CO2 = []
NO3 = []
SiO4H4 = []
totalCO2 = []
totalO2 = []

for idx,lat in enumerate(lats):
    
    global ML
    # Import the diatom only model and set its initial constraints using the Inputs.
    ML = LoadModel.LoadModel(Inputs)
    
    print("Latitude: " + str(lat))
    # Set initial values and flags. These assume that the simulation will start at daytime.
    Inputs['daycounter'] = 0
    Inputs['StartBool'] = True
    Inputs['FirstDay'] = True
    Inputs['FirstNight'] = False
    Inputs['FirstLight'] = False
    
    # Load sea-surface temperature values between March 1st and August 1st.
    # This array is at 2.5 degree latitude resolution, so ID is multiplied by 4.
    temp = np.loadtxt('.\\Data\\sstdailyMarAug2020.csv',delimiter=',')[:,idx*4]
    Henry = Extras.henry(int(temp[0]))
    CO2sat = Henry[0]*(4.21e-4)
    Inputs['Latitude'] = lat
    
    # If the solver has trouble with stiffness, the implicit Radau method can be used.
    # method='Radau',

    # Main solver. This simulation is run for two months, forcing a maximum step size
    # of one hour.
    IVPSol = scipy.integrate.solve_ivp(dydt,[0, 1440],[0.005, 0.011, CO2sat, 10, 5, 0.0005, 0, 0],max_step=1,args=[Inputs,start_date,temp,ML])
 
    time.append(IVPSol['t'])
    
    # Fill an entire array with copies of the current latitude as long as the number of timesteps.
    latinput = np.zeros((1,len(IVPSol['t'])))
    latinput.fill(lat)
    latmap.append(latinput)
    
    # Save all the simulation data.
    biomass.append(IVPSol['y'][0])
    CO2.append(IVPSol['y'][2])
    NO3.append(IVPSol['y'][3])
    SiO4H4.append(IVPSol['y'][4])
    totalBM.append(IVPSol['y'][5])
    totalCO2.append(IVPSol['y'][6])
    totalO2.append(IVPSol['y'][7])

# Data is flattened here from a 2-D array where rows represent different latitudes into a 1-D array.
timeflat = Extras.listarraytoarray(time)
latflat = Extras.flattenlat(latmap)
bioflat = Extras.listarraytoarray(biomass)
CO2flat = Extras.listarraytoarray(CO2)
NO3flat = Extras.listarraytoarray(NO3)
SiO4H4flat = Extras.listarraytoarray(SiO4H4)
totalBMflat = Extras.listarraytoarray(totalBM)
totalCO2flat = Extras.listarraytoarray(totalCO2)
totalO2flat = Extras.listarraytoarray(totalO2)

# A new 2-D array is constructed from all the flattened datasets, and saved as a .csv.
output = np.array([timeflat, latflat, bioflat, CO2flat, NO3flat, SiO4H4flat, totalBMflat, totalCO2flat, totalO2flat])
np.savetxt("surfplotallMay.csv",output,delimiter=',')

#________________________________________________
### Loading saved simulation data to a dictionary
#________________________________________________

# Re-load simulation data.
March = np.loadtxt("surfplotallMay.csv",delimiter=",")
lengths = []
lats = [55,65,75,85]

# Save the length of each solution - they will have taken a different number of time steps.
for lat in lats:
    lengths.append(len(np.where(March[1] == lat)[0]))

# Initialize arrays within a master dictionary, which will hold loaded data.
# Nine lats were previously tested, so nine rows were used.
arraydict = {}
arraydict['arraytime'] = np.zeros((9,max(lengths)))
arraydict['arraylat'] = np.zeros((9,max(lengths)))
arraydict['arraybio'] = np.zeros((9,max(lengths)))
arraydict['arrayCO2'] = np.zeros((9,max(lengths)))
arraydict['arrayNO3'] = np.zeros((9,max(lengths)))
arraydict['arraySiO4H4'] = np.zeros((9,max(lengths)))
arraydict['arraytotalBM'] =  np.zeros((9,max(lengths)))
arraydict['arraytotalCO2'] =  np.zeros((9,max(lengths)))
arraydict['arraytotalO2'] =  np.zeros((9,max(lengths)))
count = 0

# Loop through each dynamic variable:
for entry in arraydict.keys():
    
    # for each latitude, find all indexes which have the current latitude and make them their own row.
    for idx, lat in enumerate(lats):
        ids = np.where(March[1] == lat)
        hold = March[count][ids[0][0]:ids[0][-1]]
        # if this new row is shorter than the maximum row length, pad it with zeros.
        holdpad = np.pad(hold,(0,max(lengths)-len(ids[0])+1),'constant',constant_values=(0,hold[-1]))
        arraydict[entry][idx,:] = holdpad
        
    count += 1

#___________________________________________________  
### Figure 1 of Manuscript: Diatom biomass over time
#___________________________________________________

x = arraydict['arraytime']
lats = [55, 65, 75, 85]
variables = ['arraybio','arrayCO2','arrayNO3','arraySiO4H4','arraytotalBM']
colors = ['#fcaa53','#cc6a04','deepskyblue','navy']
for idy, lat in enumerate(lats):
    plt.scatter(x[idy,:],arraydict['arraybio'][idy,:],s=2,c=colors[idy])

plt.xlabel("Time (hours)")
plt.ylabel("Thalassiosira sp. Biomass (mgDW/L)")
plt.legend(["55"+chr(176)+"N","65"+chr(176)+"N","75"+chr(176)+"N","85"+chr(176)+"N"],markerscale=4)

plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=12)  # fontsize of the figure title
plt.savefig(("Figure1.pdf"),dpi=1200,format='pdf',bbox_inches='tight')

#___________________________________________________________ 
### Figure 2a of Manuscript: Nitrate concentration over time
#___________________________________________________________

x = arraydict['arraytime']
lats = ['55', '65', '75', '85']
colors = ['#fcaa53','#cc6a04','deepskyblue','navy']
for idy, lat in enumerate(lats):
    plt.scatter(x[idy,:],arraydict['arrayNO3'][idy,:],s=2,color=colors[idy],label=lat+chr(176)+"N")

plt.xlabel("Time (hours)")
plt.ylabel(u"$\mathregular{NO_3}$ Concentration (\u03bcM)") 
plt.legend(["55"+chr(176)+"N","65"+chr(176)+"N",
            "75"+chr(176)+"N","85"+chr(176)+"N"],markerscale=4)

plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=12)  # fontsize of the figure title
plt.savefig(("Figure2a.pdf"),dpi=1200,format='pdf',bbox_inches='tight')

#________________________________________________________________ 
### Figure 2b of Manuscript: Silicic acid concentration over time
#________________________________________________________________

lats = ['55', '65', '75', '85']
colors = ['#fcaa53','#cc6a04','deepskyblue','navy']
for idy, lat in enumerate(lats):
    plt.scatter(x[idy,:],arraydict['arraySiO4H4'][idy,:],s=2,color=colors[idy])

plt.xlabel("Time (hours)")
plt.ylabel(u"SiO4H4 Concentration (\u03bcM)") 
plt.ylim(bottom=0)

plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=12)  # fontsize of the figure title
plt.savefig(("Figure2b.pdf"),dpi=1200,format='pdf',bbox_inches='tight')


