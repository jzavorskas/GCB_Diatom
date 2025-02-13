"""
_____________________________________________________________
Function 1: Latitude array flattening (flattenlat) function

  Inputs: 
  
        latmap : a list of arrays with repeating latitude values,
                 each list will be the length of the number of time points at that latitude

  Outputs: 
  
        inflat : a 1-D array with all latitude values arrays stacked end to end

  To make saving data and scatter plotting easier, this function collapses latitude values into a 1-D array.
  
"""

def flattenlat(latmap):

    inflat = []
    
    for x in latmap:
        for y in list(x[0]):
            inflat.append(y)
    
    return inflat   

"""
_____________________________________________________________
Function 2: List of arrays to single array (listarraytoarray) function

  Inputs: 
  
        inarray : a list of arrays with values corresponding to separate simulations at different latitudes

  Outputs: 
  
        inflat : a 1-D array with all arrays stacked end to end

  To make scatter plotting easier, this function collapses values into a 1-D array.
  
"""

def listarraytoarray(inarray):

    inflat = []
    
    for x in inarray:
        for y in list(x):
            inflat.append(y)
    
    return inflat     

"""
_____________________________________________________________
Function 3: Adding hours for current time (addhrs) function

  Inputs: 
  
        startdate : datetime object representing the t = 0 value of simulation
        hrs : time elapsed since beginning of simulation in hours

  Outputs: 
  
        date : new datetime object representing current time

  Light intensity calculation requires date and time in the form of YYYY/MM/DD HH:MM:SS,
  this function takes the time value only in hours and generates the required form.
  
"""

def addhrs(start_date,hrs):
    from datetime import timedelta
    date = start_date + timedelta(hours=hrs)
    return date

"""
_____________________________________________________________
Function 4: Photoactive Radiation (par) calculation function

  Inputs: 
  
        DTOT : Total incident, wavelength-dependent
        allWVL : Corresponding wavelength values for DTOT

  Outputs: 
  
        par : Dictionary of lower bounds for each of the fifteen photon uptake reactions

  The par function separates the photoactive spectrum into fifteen 20 nm bins and calculates
  the molar absorption of light by Thalassiosira pseudonana by multiplying specific absorption by
  total incident radiation and converting to (umol/mg dry weight*h).
  
  This function is adapted for our model from (van Tol et al., 2021) 
"""

def par(DTOT,allWVL):
    import numpy as np
    import scipy
    # Specific absorption of light in a 20 nm bin centered at the shown wavelength.
    # i.e. "EX_photon410_e" describes the specific absorption of photons 400-420 nm.
    # UNITS: centimeter^2/mg Chlorophyll A (cm^2/mgChla)
    absorbML = {'EX_photon410_e': 4730.075289,'EX_photon430_e': 5817.128965,'EX_photon450_e':     5348.203973,'EX_photon470_e': 4050.000013,'EX_photon490_e': 3464.694801,'EX_photon510_e': 2649.794528,'EX_photon530_e': 1876.490736,'EX_photon550_e': 1334.544022,'EX_photon570_e': 873.4095179,'EX_photon590_e': 740.7816246,'EX_photon610_e': 888.7175101,'EX_photon630_e': 1082.718272,'EX_photon650_e': 1178.924274,'EX_photon670_e': 3322.974688,'EX_photon690_e': 1840.91646}
    Ipts = np.linspace(410,690,15)
    chla = 0.192e-9 # mg chlA
    DW = 16.62883198e-12 # g DW
    conv = (60*60)/1e7 #to umol mgDW-1 h-1, (60 seconds/min)(60 min/hr)(1/10000 m/cm)(1/1000 gDW/mgDW)
    
    I = []
    # Performing integration on location- and time-dependent light data to find true value of
    # incident photons in 20 nm bins from 400-700 nm.
    for idx, group in enumerate(Ipts):
        bottom = Ipts[idx] - 10
        top = Ipts[idx] + 10
        idbottom = np.flatnonzero(allWVL == bottom)[0]
        idtop = np.flatnonzero(allWVL == top)[0]
        I.append(scipy.integrate.simpson(DTOT[idbottom:idtop],allWVL[idbottom:idtop]))
    
    # Final units: umol/mgDW*hr
    pfa = {f:I[idx]*absorbML[f]*conv*(chla/DW) for idx,f in enumerate(absorbML.keys())}
    return pfa   

"""
_______________________________________________________________________________
Function 5: Temperature-dependent Henry's constant (henry) calculation function

  Inputs: 
  
        T : current temperature in degrees Celsius

  Outputs: 
  
        Henry : new datetime object representing current time

  Rubisco kinetics (below) require an accurate, temperature-dependent calculation
  of CO2 and O2 concentrations. These temperature correlations were taken from
  (Sander, 2015).
  
"""

def henry(T):
    import math
    # Conversion to Kelvin
    T += 273.15 

    # constant pressure enthalpy at reference temperature T = 298.15K. 
    # UNITS: mol/(m^3*Pa)
    H0 = [3.3e-4, # CO2
          1.3e-5, # O2
          6.4e-6]
    
    # derivative of constant pressure enthalpy versus inverse temperature
    # UNITS: Kelvin
    dlnH = [2400, # CO2
            1500, # O2
            1400] # N2
    
    Henry = [] 
    
    for i in range(len(H0)):
        
        # Units produced are mol/(m^3*Pa), converting to umol/(L*atm):
        # (101325 Pa/atm)(1,000,000 umol/mol)(1/1000 m^3/L) = 101325000
        Henryin = H0[i]*math.exp(dlnH[i]*((1/T)-(1/298.15)))*101325000
        Henry.append(Henryin)
    
    return Henry

"""
_______________________________________________________________________________
Function 6: Rubisco-limited photosynthesis of diatoms/cyanobacteria (RLPD/RLPC) calculation function

  Inputs: 
  
        CO2 : current dissolved CO2 concentration
        O2 : current dissolved O2 concentration
        T : current temperature in degrees Celsius

  Outputs: 
  
        rate : upper bound on rubisco carbon fixation rate

    Temperature- and oxygen concentration-dependent rubisco enzyme kinetics calculations.
    This function uses a competitive inhibition Michaelis-Menten expression, with temperature-dependent
    parameters (vmax, Kc, Ko). Temperature-dependent correlations were taken from (Galmes et al 2015) and 
    (Galmes et al 2016). Values at 25 degrees Celsius were taken from (Flamholz et al 2019), and conversion
    factors from kcat to vmax (Young et al 2014).
  
"""

def RLPD(CO2,O2,T):
    import math

    TK = T + 273.15 # Convert to Kelvin
    R = 8.314 # Universal gas constant, J/mol*K
    vmax25 = 10.45 # umol CO2/mgDW*h
    Ko25 = 1500 # 25degC Ko value for model diatom
    KoC3_25 = 311.6 # 25degC Ko value for C3 plant from which Ko correlation comes
    Koscale = 1500/311.6 # ratio between two Ko values
    
    ## Temperature Dependent Params for:
        ##    Kc     Ko     kcat
    scale = [ 21.1,  19.7, 30.8 ] # unitless scaling value
    dHa =   [43000,  34600, 76300 ] # J/mol
    
    
    # All temperature correlations follow the same expression:
    # parameter = e^(c-dHa/RT)    
    Kc = math.exp(scale[0] - (dHa[0]/(R*TK))) # parameters are in uM
    Ko = math.exp(scale[1] - (dHa[1]/(R*TK)))*Koscale # parameters are in uM
    kcat = math.exp(scale[2] - (dHa[2]/(R*TK)))
    vmax = kcat*vmax25 # umol/mgDW*h
    rate = (vmax*CO2)/(Kc*(1+(O2/Ko))+CO2) # umol/mgDW*h
    
    return rate

def RLPC(CO2,O2,T):
    import math
    TK = T + 273.15
    R = 8.314
    vmax25 = 3.68 # umol CO2/mgDW*h
    Ko25 = 199
    KoC3_25 = 311.6
    Koscale = 199/311.6
    ## Temperature Dependent Params for:
        ##    Kc     Ko     kcat
    scale = [ 20.8,  19.7, 14.2 ] 
    dHa =   [38800,  34600, 35200 ]
    
    ## Q10 temperature-dependence for vmax:
    
    Kc = math.exp(scale[0] - (dHa[0]/(R*TK))) # parameters are in uM
    Ko = math.exp(scale[1] - (dHa[1]/(R*TK)))*Koscale # parameters are in uM
    kcat = math.exp(scale[2] - (dHa[2]/(R*TK)))
    vmax = kcat*vmax25 # umol/mgDW*h
    rate = (vmax*CO2)/(Kc*(1+(O2/Ko))+CO2) # umol/mgDW*h
    
    return rate