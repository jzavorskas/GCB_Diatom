# Differential Equation function for IVP solver:
"""
y: [    Biomass     ] mg dry weight/L
   [    Chrysol     ] umol/L
   [ Carbon Dioxide ] umol/L
   [    Nitrate     ] umol/L
   [    Silicon     ] umol/L
   [  Total Biomass ] mg dry weight/L
   [   Total CO2    ] umol/L
   [   Total O2     ] umol/L
"""

"""
Function: Diatom-Only Ordinary Differential Equation Right-Hand-Side (DiatomOnly_dydt) calculation

  Inputs: 
  
        t : Time value in current loop iteration (hours)
        y : vector of concentrations of dynamically tracked compounds (shown above with units)
        Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations
        start_date : a datetime object representing the starting date and time, will be
                        changed throughout simulation
        temp : array with daily-resolution sea surface temperature data
        DML : COBRA model structure for diatom-only growth:
                Thalassiosira pseudonana containing all current constraints        

  Outputs: 
  
        dydt : The current rate of change of each of the dynamically tracked values,
                to be sent into an ODE solver.

This file contains the main dFBA framework, representing all calculations that occur during
a single time-step of the ODE solver. Light intensity calculation and Michaelis-Menten saturation
kinetics are calculated based on the current date and time, temperature, and nutrient concentrations.
Based on whether it is day or night, these calculations will be applied as constraints in a different way:
    Day: Photosynthesis is the only carbon source, requiring light, CO2, NO3, and Si(OH)4.
    Night: Chrysolaminarin is the only carbon source, requiring NO3, and Si(OH)4 to use and grow.

This file also diverts to other functions that solve two common dFBA problems:
    1) "fluxCalc" performs hierarchical optimization, guaranteeing a unique solution for all flux
        values used within the system of differential equations. Without this, implicit ODE solvers 
        may fail to converge.
    2) "DaySolve" and "NightSolve" perform solution saving and scaling. Often, recalculating
        the FBA linear system is unnecessary as the fluxes only change marginally, and the active
        do not change. In this situation, one can take advantage of the linear nature of FBA. 
        These two programs scale all dynamic fluxes based on the ratio between the bounds
        of the limiting nutrient (the one whose constraint is active) at the current and saved timesteps.

Learn more about these solutions here: https://doi.org/10.1002/bit.24748

"""

def DiatomOnly_dydt(t,y,Inputs,start_date,temp,DML):
    
    import math, Light, Extras, fluxCalc, DaySolve, NightSolve, scipy
    import numpy as np
    from IPython.display import clear_output 

    global daycounter
    global StartBool
    global Dfluxes
    
    # This block addresses infeasibility issues when concentrations are *very low*.
    # Slightly loosens toleranaces and institutes a quick timeout if solver is stuck.
    DML.solver.configuration.verbosity = 3
    DML.solver.configuration.timeout = 5
    DML.solver.configuration.tolerance_feasibility = 1e-6
    
    # When the date changes, date.hour will wrap back to 0. This is caught before anything else happens,
    # so all calculations should be accurate.
    elapsedday = math.floor(t/24)
    if date.hour == 0 and daycounter != elapsedday:
        daycounter += 1
        Inputs['FirstDay'] = True


    # Data for sea surface temperature is daily resolution. When the day changes, T changes.
    T = temp[elapsedday]
    
    # Resolve time; input is minutes from start, need (year, day number, HH:MM:SS)
    date = Extras.addhrs(start_date,t)
    
    # Convert date object with added time to individual values and update all light calculator values.
    Inputs['Year'] = date.year
    Inputs['DayNumber'] = date.timetuple().tm_yday
    Inputs['Hour'] = date.hour
    Inputs['Minute'] = date.minute
    Inputs['Second'] = date.second
    
    # Calculate the saturation concentration of carbon dioxide (in mM).
    Henry = Extras.henry(T)
    CO2sat = Henry[0]*(4.21e-4) # Partial pressure of CO2 in atm. # convert to uM
    O2sat = Henry[1]*(.2105) # Partial pressure of O2 in atm.
    N2sat = Henry[2]*(.781) # Partial pressure of N2 in atm.

    # Absorption spectrum of Thalassiosira pseudonana in photoactive region
    absorbML = {'EX_photon410_e': 4730.075289,'EX_photon430_e': 5817.128965,'EX_photon450_e': 5348.203973,'EX_photon470_e': 4050.000013,
            'EX_photon490_e': 3464.694801,'EX_photon510_e': 2649.794528,'EX_photon530_e': 1876.490736,'EX_photon550_e': 1334.544022,
            'EX_photon570_e': 873.4095179,'EX_photon590_e': 740.7816246,'EX_photon610_e': 888.7175101,'EX_photon630_e': 1082.718272,
            'EX_photon650_e': 1178.924274,'EX_photon670_e': 3322.974688,'EX_photon690_e': 1840.91646}

    # This block only runs the very first time-step of the ODE solver.
    # It is here to guarantee a successful first simulation regardless of
    # the external concentrations/conditions. The only constraint we set that
    # should be active is CO2 consumption.
    if StartBool:
        StartBool = False
        
        # These simulations will always start during the daytime, set startup flags accordingly.
        Inputs['NightCheck'] = False
        Inputs['FirstDay'] = True

        # Open light constraint calculation to start. No constraints applied.
        for k in absorbML.keys():
            DML.reactions.get_by_id(k).upper_bound = 0
            DML.reactions.get_by_id(k).lower_bound = -1000
        
        # Allow open fluxes for limiting nutrients to use their fluxes as flags for 
        # when the LP needs to be calculated.
        DML.reactions.EX_co2_e.lower_bound = -Extras.RLPD(CO2sat*4,O2sat,T)
        DML.reactions.EX_no3_e.lower_bound = -1000
        DML.reactions.EX_sio4h4_e.lower_bound = -1000

        # First hierarchical optimization calculation.
        [Dfluxes, Inputs, DML] = fluxCalc.fluxCalc(Inputs,DML)
        
        # Flip flags now that the first day calculation has happened.
        Inputs['FirstDay'] = False
        Inputs['FirstNight'] = True
    
    # Create an object that will store data for wavelength-dependent light transmission through the atmosphere.
    ICalc = Light.TotSolEng(Inputs)
    # Calculate the total incident extraterrestrial radiation (ETR) based on date, time and location.
    ICalc.CalcTotSolarEng(Inputs)
    # Calculate wavelength-dependent light transmission based on ETR spectrum.
    ICalc.CalcSolarSpec(Inputs)
    
    # Calculate wavelength-dependent light absorption/scattering by the ocean at a certain depth.
    ICalc.waterabs(0)
    ICalc.DTOTphoto = np.asarray(ICalc.DTOTphoto).squeeze()
    
    # Reset the night flag at every loop.
    Inputs['NightCheck'] = False
    Itot = scipy.integrate.simpson(ICalc.DTOTphoto,ICalc.allWVLphoto)
        
    # Integrate only the photoactive region (par) of the incident spectrum to get 
    # total photoactive light intensity.
    Iphoto = Extras.par(ICalc.DTOTphoto,ICalc.allWVLphoto)
    
    # Death constant under normal conditions according to (cite). 
    # Convert to h^-1.
    deathK = 1/(12*24)
    
    # If the light intensity is low, throw up the nighttime flag.
    if Itot < 10:
        Inputs['NightCheck'] = True
    
    # Running dashboard of current information from the model.
    # This is cleared every iteration of the IVP solver, so it just looks like a
    # constant dashboard.
    print("Latitude:", str(Inputs['Latitude']))
    print("Date:", date)
    print("Total Light:", Itot)
    print("Diatom Biomass:",y[0])
    print("Chrysolaminarin Concentration:",y[1])
    print("Nitrate Concentration:",y[3])
    print("Silicon Concentration:",y[4])
    print("Day Flag:", Inputs['FirstNight'])
    print("Night Flag:", Inputs['FirstDay'])
    print("Current Limiting Reaction:", Inputs['Dprevlimit'])
    print("Time Elapsed:", t)
    print(Dfluxes)
    clear_output(wait=True)

    # If the region is still in 24-hour night, this block will simply return no change
    # and assume the phytoplankton are dormant.
    if not Inputs['FirstLight']:
        if Inputs['NightCheck']:
            return [0, 0, 0, 0, 0, 0, 0, 0]
        else:
            Inputs['FirstLight'] = True
    
    # Calculate the saturation concentration of carbon dioxide (in mM).
    Henry = Extras.henry(T)
    CO2sat = Henry[0]*(4.21e-4) # Partial pressure of CO2 in atm. # convert to uM
    O2sat = Henry[1]*(.2105) # convert to uM
    
    # If the current concentration of CO2 is above its saturation concentration,
    # force it back to its saturated concentration. 
    if y[2] > CO2sat:
        y[2] = CO2sat

#___________________________________________  
### DIATOM LIMITING RESOURCES CALCULATIONS. 
#___________________________________________  
   
    # Temp-dependent Rubisco enzyme kinetics to set bounds of CO2 consumption for diatom.
    # Data can be found at:
    #   1) https://doi.org/10.1007/s11120-014-0067-8
    #   2) https://doi.org/10.1093/jxb/erw267
    DML.reactions.EX_co2_e.upper_bound = 1000
    DML.reactions.EX_co2_e.lower_bound = -Extras.RLPD(CO2sat*4,O2sat,T)

    # Conversion factors necessary for certain calculations involving biomass.
    # Data can be found at: https://doi.org/10.1371/journal.pone.0241960
    chla = 0.192e-9
    DW = 16.62883198e-12

#    ML.reactions.EX_hco3_e.lower_bound = -((2.438*hco3)/(258.6+hco3)) # Trimborn, 2009 [Not Used]
    DML.reactions.EX_hco3_e.lower_bound = 0

    # Michaelis-Menten saturation kinetics for nitrate uptake. To avoid errors,
    # NO3 concentration is set to zero when it becomes severely limiting.
    # Low but not zero NO3 concentrations can cause infeasible results.
    # Data can be found at: https://doi.org/10.3354/meps09088
    if y[3] > 0.01:
        DML.reactions.EX_no3_e.lower_bound = -((1.474*y[3])/(6.14+y[3])) 
    else:
        y[3] = 0
        DML.reactions.EX_no3_e.lower_bound = 0

    # Modified Michaelis-Menten saturation kinetics for silicic acid uptake, using a Hill coefficient
    # which represents cooperative allosteric interaction during silicic acid uptake. 
    # To avoid errors, Si(OH)4 concentration is set to zero when it becomes severely limiting.
    # Low but not zero Si(OH)4 concentrations can cause infeasible results.
    # Data can be found at: https://doi.org/10.1104/pp.107.107094
    n = 1.9
    if y[4] > 0.1:   
        DML.reactions.EX_sio4h4_e.lower_bound = -((1.961*(y[4]**n))/((8.1**n)+(y[4]**n)))
    else:    
        y[4] = 0
        DML.reactions.EX_sio4h4_e.lower_bound = 0

    # Defaults for light intensity are ub = 0, lb = 0. 
    
    # If it is nighttime, carbon fixation is set to zero and only stored carbon
    # can be used for respiration.
    if Inputs['NightCheck']: 
        # Stored carbon, in the form of chrysolaminarin is consumed at night.
        # See supplemental files for our derivation of a kinetic model for chrysolaminarin consumption. 
        # Original kinetic data is from: https://doi.org/10.1016/j.ijbiomac.2023.126361
        DML.reactions.EX_chryso_e.lower_bound = -(0.935*y[1])/(133.24+y[1]) 
        # No light, no photosynthesis.
        DML.reactions.EX_co2_e.lower_bound = 0

        # Allow a very small amount of light uptake. Not enough for photosynthesis to occur, but
        # enough to allow creation of necessary micronutrients for stored carbon metabolism.
        for k in Iphoto.keys():
            DML.reactions.get_by_id(k).upper_bound = 0
            DML.reactions.get_by_id(k).lower_bound = -1
            DML.reactions.get_by_id(k).upper_bound = -0.9999   

    # If it is daytime, we assume stored carbon is not consumed at all, and fixed carbon
    # is the sole food source.        
    else:
        DML.reactions.EX_chryso_e.lower_bound = 0
        for k in Iphoto.keys():
            DML.reactions.get_by_id(k).upper_bound = 0
            DML.reactions.get_by_id(k).lower_bound = Iphoto[k]*-1.
            DML.reactions.get_by_id(k).upper_bound = Iphoto[k]*-0.9999

#_______________________________________________   
### LP Calculation Flags:
# abort loop here if no flags are triggered and use previous LP results.
# flags typically only trigger when a day/night change occurs.
#_______________________________________________
    
    
    # The first time night occurs, a basis change is necessary.
    # Calculate fluxes and flip day/night flags.
    if Inputs['NightCheck'] and Inputs['FirstNight']:
        print("Basis Change at ", t)
        [Dfluxes, Inputs, DML] = fluxCalc.fluxCalc(Inputs,DML)
        Inputs['Dprevflux'] = Dfluxes
        Inputs['FirstNight'] = False
        Inputs['FirstDay'] = True    

    # Any other time night occurs, recalculation is likely not necessary!
    elif Inputs['NightCheck'] and not Inputs['FirstNight']:
        [Dfluxes, Inputs, DML] = NightSolve.NightSolve(t,y,Inputs,Dfluxes,DML)

    # The first time day occurs, a basis change is necessary.
    # Calculate fluxes and flip day/night flags.
    elif not Inputs['NightCheck'] and Inputs['FirstDay']:
        print("Basis Change at ", t)
        [Dfluxes, Inputs, DML] = fluxCalc.fluxCalc(Inputs,DML)
        Inputs['prevflux'] = Dfluxes
        Inputs['FirstDay'] = False
        Inputs['FirstNight'] = True

    # Any other time day occurs, recalculation is likely not necessary!    
    elif not Inputs['NightCheck'] and not Inputs['FirstDay']:
        [Dfluxes, Inputs, DML] = DaySolve.DaySolve(t,y,Inputs,Dfluxes,DML)
    elif not Inputs['NightCheck'] and Inputs['Dprevlimit'] == "EX_chryso_e":
        [Dfluxes, Inputs, DML] = fluxCalc.fluxCalc(Inputs,DML)        
    
#_______________________________________________   
### ODE Updating Section:
#_______________________________________________

    # Check for nighttime!
    if Inputs['NightCheck']:
        
        if y[1] > 0 and y[3] > 1e-2 and y[4] > 1e-1:
            
            dydt = [(Dfluxes["DM_biomass_c"]-deathK)*y[0]] # Biomass
            dydt.append(Dfluxes["EX_chryso_e"]*y[0]) # Chyrsolaminarin
            dydt.append(0) # Carbon Dioxide
            dydt.append(Dfluxes["EX_no3_e"]*y[0])
            dydt.append(Dfluxes["EX_sio4h4_e"]*y[0])            
            dydt.append(Dfluxes["DM_biomass_c"]*y[0]) # Biomass without death constant
            dydt.append(Dfluxes["EX_o2_e"]*y[0])
            dydt.append(0) # Formate respiration acceptor, no O2 produced
        
        # When concentrations are too low, diatoms slowly die off.
        else:
            dydt = [-deathK*y[0], 0, 0, 0, 0, 0, 0, 0]
            return dydt
            
    # Otherwise, it is daytime.
    else:
                
        dydt = [(Dfluxes["DM_biomass_c"]-deathK)*y[0]] # Live biomass
        if (Dfluxes["DM_biomass_c"]-deathK)*y[0] > 0:
            dydt.append((0.0085*math.exp(0.0163*Itot))*y[0])   # Chrysolaminarin
        else:
            dydt.append(0)
        dydt.append(0) # Carbon Dioxide
        dydt.append(Dfluxes["EX_no3_e"]*y[0]) # No nitrate diffusion if they are dead.
        dydt.append(Dfluxes["EX_sio4h4_e"]*y[0])  
        dydt.append(Dfluxes["DM_biomass_c"]*y[0]) # Biomass without death constant
        dydt.append(Dfluxes["EX_co2_e"]*y[0])
        dydt.append(Dfluxes["EX_o2_e"]*y[0])
    
    return dydt