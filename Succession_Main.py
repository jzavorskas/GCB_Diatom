"""
Author: Joe Zavorskas
Date Started: 6/26/2023
Last Edit: 04/11/2024

Main file for diatom succession simulations. Simulates the annual succession pattern between
Thalassiosira species and Chaetoceros species in the Arctic. Symbiosis between Chaetoceros and 
a cyanobacterium is included, which relieves nitrogen stress on Chaetoceros after nitrate has been depleted.

A detailed input instance is passed throughout the differential equation function defined below, which 
contains information required to calculate light intensity, and saved LP solutions for the diatom-only FBA
and the symbiotic FBA problems.

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
from Succession_dydt import Succession_dydt as dydt

# Initialize the inputs dictionary which will be used through dFBA simulation.
[Inputs, start_date] = Input.Input()

#__________________________
### Main Simulation Section
#__________________________

# Initialize all arrays that will hold dFBA data.
lats = [55, 65, 75, 85]
time = []
latmap = []
TPbiomass = []
CSbiomass = []
Cybiomass = []
TPchryso = []
CSchryso = []
NO3 = []
SiO4H4 = []
Fe = []
totalBM = []
totalCO2 = []
totalO2 = []
totalN2 = []

for idx,lat in enumerate(lats):
    
    global ML
    # Import three models, and set their constraints to default:
    # 1) diatom only, 2) symbiotic, 3) cyanobacteria only, for reference.
    [ML, DML, refmodel, IYSi, INSi, Inputs] = MergeLoadModel.MergeLoadModel(Inputs)
    
    # Store the interpolators inside of the Inputs dictionary.
    Inputs['IYSi'] = IYSi
    Inputs['INSi'] = INSi

    print("Latitude: " + str(lat))
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
    
    # If the solver has trouble with stiffness, the implicit Radau method can be used.
    # method='Radau',

    # Main solver. This simulation is run for two months, forcing a maximum step size
    # of one hour.
    IVPSol = scipy.integrate.solve_ivp(dydt,[0, 3600],[0.005, 0.005, 0.0008, 0.011, 0.0011, 9, 4, 0.004, 0.0105, 0, 0, 0],max_step=1,args=[Inputs,start_date,temp,ML,DML,refmodel])
 
    time.append(IVPSol['t'])
    
    # Fill an entire array with copies of the current latitude as long as the number of timesteps.
    latinput = np.zeros((1,len(IVPSol['t'])))
    latinput.fill(lat)
    latmap.append(latinput)

    # Save all the simulation data.
    TPbiomass.append(IVPSol['y'][0])
    CSbiomass.append(IVPSol['y'][1])
    Cybiomass.append(IVPSol['y'][2])
    TPchryso.append(IVPSol['y'][3])
    CSchryso.append(IVPSol['y'][4])
    NO3.append(IVPSol['y'][5])
    SiO4H4.append(IVPSol['y'][6])
    Fe.append(IVPSol['y'][7])
    totalBM.append(IVPSol['y'][8])
    totalCO2.append(IVPSol['y'][9])
    totalO2.append(IVPSol['y'][10])
    totalN2.append(IVPSol['y'][11])

# Data is flattened here from a 2-D array where rows represent different latitudes into a 1-D array.    
timeflat = Extras.listarraytoarray(time)
latflat = Extras.flattenlat(latmap)
TPbioflat = Extras.listarraytoarray(TPbiomass)
CSbioflat = Extras.listarraytoarray(CSbiomass)
Cybioflat = Extras.listarraytoarray(Cybiomass)
TPchrysoflat = Extras.listarraytoarray(TPchryso)
CSchrysoflat = Extras.listarraytoarray(CSchryso)
NO3flat = Extras.listarraytoarray(NO3)
SiO4H4flat = Extras.listarraytoarray(SiO4H4)
Feflat = Extras.listarraytoarray(Fe)
totalBMflat = Extras.listarraytoarray(totalBM)
totalCO2flat = Extras.listarraytoarray(totalCO2)
totalO2flat = Extras.listarraytoarray(totalO2)
totalN2flat = Extras.listarraytoarray(totalN2)

# A new 2-D array is constructed from all the flattened datasets, and saved as a .csv.
output = np.array([timeflat, latflat, TPbioflat, CSbioflat, Cybioflat, TPchrysoflat, CSchrysoflat, 
                   NO3flat, SiO4H4flat, Feflat, totalBMflat, totalCO2flat, totalO2flat, totalN2flat])
np.savetxt("SuccessionMarch2.csv",output,delimiter=',')

#________________________________________________
### Loading saved simulation data to a dictionary
#________________________________________________

# Re-load simulation data.
March = np.loadtxt("SuccessionMarch2.csv",delimiter=",")
lengths = []
lats = [55, 65, 75, 85]

# Save the length of each solution - they will have taken a different number of time steps.
for lat in lats:
    lengths.append(len(np.where(March[1] == lat)[0]))

# Initialize arrays within a master dictionary, which will hold loaded data.
# Nine lats were previously tested, so nine rows were used.
arraydict = {}
arraydict['arraytime'] = np.zeros((9,max(lengths)))
arraydict['arraylat'] = np.zeros((9,max(lengths)))
arraydict['arrayTPbio'] = np.zeros((9,max(lengths)))
arraydict['arrayCSbio'] = np.zeros((9,max(lengths)))
arraydict['arrayCybio'] = np.zeros((9,max(lengths)))
arraydict['arrayTPchryso'] = np.zeros((9,max(lengths)))
arraydict['arrayCSchryso'] = np.zeros((9,max(lengths)))
arraydict['arrayNO3'] = np.zeros((9,max(lengths)))
arraydict['arraySiO4H4'] = np.zeros((9,max(lengths)))
arraydict['arrayFe'] = np.zeros((9,max(lengths)))
arraydict['arraytotalBM'] =  np.zeros((9,max(lengths)))
arraydict['arraytotalCO2'] =  np.zeros((9,max(lengths)))
arraydict['arraytotalO2'] =  np.zeros((9,max(lengths)))
arraydict['arraytotalN2'] =  np.zeros((9,max(lengths)))
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
colors = ['#fcaa53','#cc6a04','deepskyblue','navy']
for idy, lat in enumerate(lats):
    plt.scatter(x[idy,:],arraydict['arrayTPbio'][idy,:],s=0.75,c=colors[idy])

plt.xlabel("Time (hours)")
plt.ylabel("Thalassiosira sp. Biomass (mgDW/L)")
plt.legend(["55"+chr(176)+"N","65"+chr(176)+"N","75"+chr(176)+"N","85"+chr(176)+"N"],markerscale=6)

plt.rc('font', size=15)          # controls default text sizes
plt.rc('axes', titlesize=15)     # fontsize of the axes title
plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
plt.rc('legend', fontsize=8)    # legend fontsize
plt.rc('figure', titlesize=12)  # fontsize of the figure title
plt.savefig(("Figure1.pdf"),dpi=1200,format='pdf')

#____________________________________________________________  
### Figure 5 of Manuscript: Combined diatom biomass over time
#____________________________________________________________

x = arraydict['arraytime']
colors = ['#fcaa53','#cc6a04','deepskyblue','navy']
for idy, lat in enumerate(lats):
    plt.scatter(x[idy,:],(arraydict['arrayTPbio'][idy,:]+arraydict['arrayCSbio'][idy,:]),s=0.75,c=colors[idy])

plt.xlabel("Time (hours)")
plt.ylabel("Combined Diatom Biomass (mgDW/L)")
plt.legend(["55"+chr(176)+"N","65"+chr(176)+"N","75"+chr(176)+"N","85"+chr(176)+"N"],markerscale=6)

plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=12)  # fontsize of the figure title
plt.savefig(("SymFigure5.pdf"),dpi=1200,format='pdf',bbox_inches='tight')

#_____________________________________________________________________  
### Figure 4 of Manuscript: Symbiotic diatom biomass over time ONLY
#_____________________________________________________________________

x = arraydict['arraytime']
colors = ['#fcaa53','#cc6a04','deepskyblue','navy']
for idy, lat in enumerate(lats):
    plt.scatter(x[idy,:],arraydict['arrayCSbio'][idy,:],s=0.75,c=colors[idy])

plt.xlabel("Time (hours)")
plt.ylabel("Chaetoceros sp. Biomass (mgDW/L)")
plt.legend(["55"+chr(176)+"N","65"+chr(176)+"N","75"+chr(176)+"N","85"+chr(176)+"N"],markerscale=8)

plt.rc('font', size=15)          # controls default text sizes
plt.rc('axes', titlesize=15)     # fontsize of the axes title
plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
plt.rc('legend', fontsize=20)    # legend fontsize
plt.rc('figure', titlesize=12)  # fontsize of the figure title
plt.savefig(("SymFigure4.pdf"),dpi=1200,format='pdf',bbox_inches='tight')

#______________________________________________________________________________________  
### Not a figure, just for fun: separated diatom biomass over time, to show succession
#______________________________________________________________________________________

x = arraydict['arraytime']
lats = [55, 65, 75, 85]
variables = ['arrayTPbio','arrayCSbio','arrayCybio','arrayNO3','arraySiO4H4','arrayFe','arraytotalBM','arraytotalCO2','arraytotalO2','arraytotalN2']
labels = ["Active T.pseudonana Biomass \n (mgDW/L)", "Active C.socialis Biomass \n (mgDW/L)", "Active Cyanobacteria Biomass \n (mgDW/L)", "Nitrate Concentration \n (umol/L)",
           "Silicic Acid Concentration \n (umol/L)", "Iron Concentration \n (nmol/L)", "Total Biomass Generated \n (mgDW/L)", "Total CO2 Consumed \n (umol/L)",
           "Total O2 Produced \n (umol/L)", "Total N2 Consumed \n (umol/L)"]
filenames = ["TPBio","CSBio","CyanoBio","Nitrate","Silicic","Iron","TotalBio","TotalCO2","TotalO2","TotalN2"]

for idx, name in enumerate(filenames):
    plt.figure()

for idy, lat in enumerate(lats):
    plt.scatter(x[idy,:],arraydict['arrayTPbio'][idy,:],s=0.5)
    plt.scatter(x[idy,:],arraydict['arrayCSbio'][idy,:],s=0.5)

plt.xlabel("Time (hours)")
plt.ylabel(labels[idx]) 
plt.legend(["55"+chr(176)+"N TP","55"+chr(176)+"N CS","65"+chr(176)+"N TP","65"+chr(176)+"N CS",
            "75"+chr(176)+"N TP","75"+chr(176)+"N CS","85"+chr(176)+"N TP","85"+chr(176)+"N CS"],markerscale=6)

plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=8)    # legend fontsize
plt.rc('figure', titlesize=12)  # fontsize of the figure title
plt.savefig(("SuccessionSeparateBio.pdf"),dpi=1200,format='pdf')

#________________________________________________________________________________________  
### Figure 6 of Manuscript: Limiting nutrients (nitrogen, silicic acid, iron) over time
#________________________________________________________________________________________

colors = ['blue','green','red','purple']
latstr = ['55', '65', '75', '85']
FigureMarker = ['a','b','c','d']
for idy, lat in enumerate(lats):
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Time (hours)')
    ax1.set_ylabel(u"$\mathregular{NO_3}$/$\mathregular{Si(OH)_4}$ Concentration (\u03bcM)")
    ax1.plot(x[idy,:2500],arraydict['arrayNO3'][idy,:2500],color='gold',label="Nitrate")    
    ax1.plot(x[idy,:2500],arraydict['arraySiO4H4'][idy,:2500],color='navy',label="Silicon")
    ax1.set_ylim(bottom=0)
    # ax1.legend(loc="upper right")
    # plt.scatter(x[idy,:],arraydict['arrayNO3'][idy,:],color=colors[idy],label=latstr[idy]+chr(176)+"N",s=0.75)    
    # plt.scatter(x[idy,:],arraydict['arraySiO4H4'][idy,:],color=colors[idy],s=0.75)
    # plt.scatter(x[idy,:],arraydict['arrayFe'][idy,:],color=colors[idy],s=0.75)
    
    ax2 = ax1.twinx()
    
    ax2.set_ylabel('Iron Concentration (nM)')
    ax2.plot(x[idy,:2500],arraydict['arrayFe'][idy,:2500]*1000,color='#cc6a04',label="Iron")
    # ax2.legend(loc="upper right",bbox_to_anchor=(1,0.86))
    
    fig.legend(["$\mathregular{NO_3}$","$\mathregular{Si(OH)_4}$","Fe"],markerscale=10, bbox_to_anchor=(0.9,0.88))
    
    # plt.xlabel("Time (hours)")
    # plt.ylabel("NO3/SiO4H4 Concentration \n (umol)") 
    # plt.legend(["Nitrate","Silicon","Iron"],markerscale=8)

    plt.rc('font', size=16)          # controls default text sizes
    plt.rc('axes', titlesize=14)     # fontsize of the axes title
    plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
    plt.rc('legend', fontsize=12)    # legend fontsize
    plt.rc('figure', titlesize=12)  # fontsize of the figure title
    # fig.tight_layout()
    plt.savefig(("SymFigure6"+FigureMarker[idy]+".pdf"),dpi=1200,format='pdf',bbox_inches='tight')