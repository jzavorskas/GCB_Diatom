"""
Description Required
"""

# Import necessary packages:

import cobra
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
from datetime import datetime, timedelta
from Light import TotSolEng
from scipy.interpolate import RegularGridInterpolator
import netCDF4
from IPython.display import clear_output

import MergeDaySolve, MergeNightSolve, MergefluxCalc, MergeLoadModel, Light, Extras

### Section Break Here

# Initialize a dictionary to hold all inputs to the class.
Inputs = {}

# User defined inputs here!!! Names will be more descriptive in this user-facing section.
Inputs['Year'] = 2023
Inputs['DayNumber'] = 0
Inputs['Month'] = 3 # Month 0 is a flag for not using MM/DD/YYYY mode.
Inputs['Day'] = 1 # Day 0 is a flag for not using MM/DD/YYYY mode.
Inputs['Hour'] = 12 # 24-hour time
Inputs['Minute'] = 0
Inputs['Second'] = 0
Inputs['Interval'] = 0 # Defaults to zero. Do not change, interval will be handed in dFBA.
Inputs['Latitude'] = 75 # Bottom of Arctic Circle
Inputs['Longitude'] = 0 # Prime Meridian
Inputs['Timezone'] = 0 # We are at the prime meridian. (This will likely be automatically calculated.)
Inputs['Pressure'] = 1013 # DEFAULT: 1013 mb = 1 atmosphere
Inputs['Temperature'] = 0 # DEFAULT: 10 deg C (arctic will be much colder)

Inputs['Turbidity'] = 0.084
Inputs['Water Vapor'] = 1.4164
Inputs['Ozone'] = 0.3438
Inputs['Albedo'] = 0.1

# Input wavelength-dependent data tables.
Inputs['WVL-ETR'] = np.loadtxt('./Data/WVL-ETR.csv',delimiter=',',dtype='float64', encoding='utf-8-sig')
Inputs['WVL-ABS'] = np.loadtxt('./Data/WVL-ABS.csv',delimiter=',',dtype='float64', encoding='utf-8-sig')

Inputs["NO3flag"] = 0
Inputs["SiOH4flag"] = 0

start = str(Inputs['Month']).zfill(2) + "/" + str(Inputs['Day']).zfill(2) + "/" + str(Inputs['Year']) + " " + str(Inputs['Hour'])
start_date = datetime.strptime(start, "%m/%d/%Y %H")

### Section Break Here

global daycounter
global StartBool
global FirstDay
global FirstNight
global FirstLight

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

    ML, refmodel, IYSi, INSi = LoadModel()
    
    print("Latitude: " + str(lat))
    daycounter = 0
    StartBool = True
    FirstDay = True
    FirstNight = False
    FirstLight = False
    
    temp = np.loadtxt('sstdailyMarAug2020.csv',delimiter=',')[:,idx*4]
    Henry = henry(int(temp[0]))
    CO2sat = Henry[0]*(4.21e-4)
    Inputs['Latitude'] = lat

    refmodel.reactions.get_by_id("EX_cpd00011[e]").lower_bound = -RLPC(CO2sat*4,O2sat,5)
    refmodel.reactions.get_by_id("EX_cpd00242[e]").lower_bound = 0
    refmodel.reactions.get_by_id("EX_cpd00027[e]").lower_bound = 0
    refmodel.reactions.biomass_eq_33047__hc.lower_bound = 0
    
    IVPSol = scipy.integrate.solve_ivp(dydt,[0, 4300],[0.0005, 0.0005, 0.055, 27, 10, 10, 0.001, 0, 0, 0],max_step=1,args=[Inputs,start_date,temp,IYSi,INSi,Tn,datamin])
    #IVPSol = scipy.integrate.solve_ivp(dydt,[0, 720],[7.472547668746805, 8.554385159333298, 180.5802403287163, 27, 0, 0.28245253147206567, 0, 0, 0, 0],max_step=1,args=[Inputs,start_date,temp,IAll,Tn,datamin])
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
    
timeflat = listarraytoarray(time)
latflat = flattenlat(latmap)
Dbioflat = listarraytoarray(Dbiomass)
Cbioflat = listarraytoarray(Cbiomass)
NO3flat = listarraytoarray(NO3)
SiO4H4flat = listarraytoarray(SiO4H4)
totalBMflat = listarraytoarray(totalBM)
totalCO2flat = listarraytoarray(totalCO2)
totalO2flat = listarraytoarray(totalO2)
totalN2flat = listarraytoarray(totalN2)

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



