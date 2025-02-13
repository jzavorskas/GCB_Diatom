# Differential Equation function for IVP solver:
"""
y: 
0  [    T.pseudo Biomass   ] mg dry weight/L
1  [   C.socialis Biomass  ] mg dry weight/L
2  [    Cyano Biomass      ] mg dry weight/L
3  [    T.pseudo Chryso    ] umol/L
4  [   C.socialis Chryso   ] umol/L
5  [       Nitrate         ] umol/L
6  [       Silicon         ] umol/L
7  [        Iron           ] nmol/L
8  [    Total Biomass      ] mg dry weight/L
9  [      Total CO2        ] umol/L
10 [      Total O2         ] umol/L
11 [      Total N2         ] umol/L
"""

"""
Function: Diatom Succession Ordinary Differential Equation Right-Hand-Side (Succession_dydt) calculation

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
                represents Thalassiosira-type diatoms that dominate first in succession
        ML :  COBRA model structure for symbiotic growth:
                Thalassiosira pseudonana containing all current constraints
                represents Chaetoceros-type diatoms that take over in nitrogen starvation   
        refmodel : COBRA model structure representing cyanobacteria-only growth
                    Anabaena variabilis containing all current constraints
                    This model is used to find basal cyanobacteria metabolism 
                    (what it would do if not in symbiosis) 

  Outputs: 
  
        dydt : The current rate of change of each of the dynamically tracked values,
                to be sent into an ODE solver.

This file contains the main dFBA framework, representing all calculations that occur during
a single time-step of the ODE solver. Light intensity calculation and Michaelis-Menten saturation
kinetics are calculated based on the current date and time, temperature, and nutrient concentrations.
Based on whether it is day or night, these calculations will be applied as constraints in a different way:
    Day: Photosynthesis is the only carbon source, requiring light, CO2, NO3, Si(OH)4, and Fe.
    Night: Chrysolaminarin is the only carbon source, requiring NO3, and Si(OH)4 to metabolize carbon and grow.

The diatom-only model (DML) and symbiotic model (ML) directly compete for resources, pulling
from a common pool of nitrogen, silicon, and iron.

This file also diverts to other functions that solve two common dFBA problems:
    1) "MergefluxCalc" and "fluxCalc" perform hierarchical optimization, guaranteeing a unique solution for all flux
        values used within the system of differential equations. Without this, implicit ODE solvers 
        may fail to converge. These two functions are used separately to optimize diatom-only and symbiotic, respectively.
    2) "MergeDaySolve/DaySolve" and "MergeNightSolve/NightSolve" perform solution saving and scaling. 
        Often, recalculating the FBA linear system is unnecessary as the fluxes only change marginally, 
        and the active constraints do not change. In this situation, one can take advantage of the linear nature of FBA. 
        These two programs scale all dynamic fluxes based on the ratio between the bounds
        of the limiting nutrient (the one whose constraint is active) at the current and saved timesteps.
        The "Merge" version of each file simply adds in more dynamic fluxes that must be tracked (i.e. iron).

Learn more about these solutions here: https://doi.org/10.1002/bit.24748

"""

def Succession_dydt(t,y,Inputs,start_date,temp,ML,DML,refmodel):
    
    import math, scipy
    from IPython.display import clear_output
    import numpy as np
    from Code_NoDepth.Hierarchical_Optimization import MergefluxCalc, MergeDaySolve, MergeNightSolve
    from Code_NoDepth.Hierarchical_Optimization import fluxCalc, DaySolve, NightSolve
    from Code_NoDepth.Auxiliary_Functions import Light, Extras, Input, LoadModel, MergeLoadModel
    # Auxiliary Functions and Hierarchical Optimization MUST be in the same parent directory
    # as this file for import to work.

    global fluxes
    global Dfluxes
    
    # This block addresses infeasibility issues when concentrations are *very low*.
    # Slightly loosens toleranaces and institutes a quick timeout if solver is stuck.
    # ML.solver.configuration.verbosity = 3
    ML.solver.configuration.timeout = 5
    ML.solver.configuration.tolerance_feasibility = 1e-6
    
    # Resolve time; input is minutes from start, need (year, day number, HH:MM:SS)
    date = Extras.addhrs(start_date,t)

    # When the date changes, date.hour will wrap back to 0. This is caught before anything else happens,
    # so all calculations should be accurate.
    elapsedday = math.floor(t/24)

    # Data for sea surface temperature is daily resolution. When the day changes, T changes.
    T = temp[elapsedday]

    # Calculate the saturation concentration of carbon dioxide (in mM).
    Henry = Extras.henry(T)
    CO2sat = Henry[0]*(4.21e-4) # Partial pressure of CO2 in atm.
    O2sat = Henry[1]*(.2105) # Partial pressure of O2 in atm.
    N2sat = Henry[2]*(.781) # Partial pressure of N2 in atm.
    Inputs['sat'] = [CO2sat,O2sat,N2sat]

    # Sea surface temperature data is daily. At midnight, update the reference model for cyanobacteria with
    # the new temperature data.
    if date.hour == 0 and Inputs['daycounter'] != elapsedday:
        refmodel.reactions.get_by_id("EX_cpd00011[e]").lower_bound = -Extras.RLPC(CO2sat*4,O2sat,T)
        Inputs['refsolution'] = refmodel.optimize()
        Inputs['LightLimit'] = Inputs['refsolution'].fluxes["EX_PHO1"]+Inputs['refsolution'].fluxes["EX_PHO2"]
        Inputs['NH4Constraint'] = -4.5*Inputs['refsolution'].fluxes["EX_cpd00528[e]"]
        
        Inputs['daycounter'] += 1
        Inputs['FirstDay'] = True

    if date.hour == 0 and Inputs['daycounter'] != elapsedday:
        Inputs['daycounter'] += 1
        Inputs['FirstDay'] = True
    
    # Convert date object with added time to individual values and update all light calculator values.
    Inputs['Year'] = date.year
    Inputs['DayNumber'] = date.timetuple().tm_yday
    Inputs['Hour'] = date.hour
    Inputs['Minute'] = date.minute
    Inputs['Second'] = date.second
    
    # Absorption spectrum of Thalassiosira pseudonana in photoactive region
    absorbML = {'EX_photon410_e': 4730.075289,'EX_photon430_e': 5817.128965,'EX_photon450_e': 5348.203973,'EX_photon470_e': 4050.000013,
            'EX_photon490_e': 3464.694801,'EX_photon510_e': 2649.794528,'EX_photon530_e': 1876.490736,'EX_photon550_e': 1334.544022,
            'EX_photon570_e': 873.4095179,'EX_photon590_e': 740.7816246,'EX_photon610_e': 888.7175101,'EX_photon630_e': 1082.718272,
            'EX_photon650_e': 1178.924274,'EX_photon670_e': 3322.974688,'EX_photon690_e': 1840.91646}
    
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
        
    # Integrate only the photoactive region of the incident spectrum to get total photoactive light intensity.
    Iphoto = Extras.par(ICalc.DTOTphoto,ICalc.allWVLphoto)
    
    # Death constant under normal conditions according to (cite). 
    # Convert to h^-1.
    deathK = 1/(12*24)
    
    # If the light intensity is low, throw up the nighttime flag.
    if Itot < -2*Inputs["LightLimit"]:
        Inputs['NightCheck'] = True
    
    # Running dashboard of current information from the model.
    # This is cleared every iteration of the IVP solver, so it just looks like a
    # constant dashboard.
    if not Inputs['StartBool']:
        print("Latitude:", str(Inputs['Latitude']))
        print("Date:", date)
        print("Total Light:", Itot)
        print("T. pseudonana Biomass:",y[0])
        print("C. socialis Biomass:",y[1])
        print("A. variabilis Biomass:",y[2])
        print("T.pseudo Chrysolaminarin:",y[3])
        print("C.social Chrysolaminarin:",y[3])
        print("Nitrate Concentration:",y[5])
        print("Silicon Concentration:",y[6])
        print("Iron Concentration:",y[7])
        print("Day Flag:", Inputs['FirstNight'])
        print("Night Flag:", Inputs['FirstDay'])
        print("Current Limiting Reaction:", Inputs['prevlimit'])
        print("Time Elapsed:", t)
        print(Dfluxes)
        print(fluxes)
        clear_output(wait=True)

    # If the region is still in 24-hour night, this block will simply return no change
    # and assume the phytoplankton are dormant.
    if not Inputs['FirstLight']:
        if Itot < -2*Inputs["LightLimit"]:
            return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        else:
            Inputs['FirstLight'] = True
    
    # Calculate the saturation concentration of carbon dioxide (in mM).
    Henry = Extras.henry(T)
    CO2sat = Henry[0]*(4.21e-4) # Partial pressure of CO2 in atm. # convert to uM
    O2sat = Henry[1]*(.2105) # convert to uM
    N2sat = Henry[2]*(.781)

    # If the current concentration of CO2 is above its saturation concentration,
    # force it back to its saturated concentration. 
    if y[3] > CO2sat:
        y[3] = CO2sat

    # "Upper Bound" - CO2 production is left open, typically doesn't happen.
    # "Lower Bound" - RuBisCo-limited CO2 uptake kinetics
    # Temp-dependent Rubisco enzyme kinetics to set bounds of CO2 consumption for diatom.
    # Data can be found at:
    #   1) https://doi.org/10.1007/s11120-014-0067-8
    #   2) https://doi.org/10.1093/jxb/erw267
    ML.reactions.EX_co2_e.upper_bound = 1000
    ML.reactions.EX_co2_e.lower_bound = -Extras.RLPD(CO2sat*4,O2sat,T)
    DML.reactions.EX_co2_e.upper_bound = 1000
    DML.reactions.EX_co2_e.lower_bound = -Extras.RLPD(CO2sat*4,O2sat,T)
    
    # Bicarbonate uptake kinetics are turned off. RuBisCo can only consume CO2, but HCO3 is converted
    # to CO2 by carbonic anhydride in diatom (and cyanobacteria) within their carbon-concentrating
    # mechanisms. This effect is approximated by multiplying CO2's saturating concentration by 4, a
    # conservative estimate for pyrenoid CO2 concentration from:
    # www.pnas.org/cgi/doi/10.1073/pnas.1018062108
    ML.reactions.EX_hco3_e.lower_bound = 0
    DML.reactions.EX_hco3_e.lower_bound = 0
    
    # Michaelis-Menten saturation kinetics for nitrate uptake. To avoid errors,
    # NO3 concentration is set to zero when it becomes severely limiting.
    # Low but not zero NO3 concentrations can cause infeasible results.
    # Data can be found at: https://doi.org/10.3354/meps09088
    if y[5] > 0.01:
        ML.reactions.EX_no3_e.lower_bound = -((1.474*y[5])/(6.14+y[5]))
        DML.reactions.EX_no3_e.lower_bound = -((1.474*y[5])/(6.14+y[5]))
    else:
        y[5] = 0
        ML.reactions.EX_no3_e.lower_bound = 0
        DML.reactions.EX_no3_e.lower_bound = 0
    
    # Modified Michaelis-Menten saturation kinetics for silicic acid uptake, using a Hill coefficient
    # which represents cooperative allosteric interaction during silicic acid uptake. 
    # To avoid errors, Si(OH)4 concentration is set to zero when it becomes severely limiting.
    # Low but not zero Si(OH)4 concentrations can cause infeasible results.
    # Data can be found at: https://doi.org/10.1104/pp.107.107094
    n = 1.9
    if y[6] > 0.1:
        ML.reactions.EX_sio4h4_e.lower_bound = -((1.961*(y[6]**n))/((8.1**n)+(y[6]**n)))
        DML.reactions.EX_sio4h4_e.lower_bound = -((1.961*(y[6]**n))/((8.1**n)+(y[6]**n)))
    else:    
        y[6] = 0
        ML.reactions.EX_sio4h4_e.lower_bound = 0
        DML.reactions.EX_sio4h4_e.lower_bound = 0

    # Michaelis-Menten saturation kinetics for siderophore-mediated iron uptake
    # by cyanobacteria.
    # To avoid errors, Fe concentration is set to zero when it becomes severely limiting.
    # Low but not zero Fe concentrations can cause infeasible results.
    # Data can be found at: https://doi.org/10.1111/mmi.13049
    if y[7] > 5e-5:
        ML.reactions.get_by_id("EX_cpd10515[e]").lower_bound = -((0.0024*y[7])/(0.0023+y[7]))
    else:
        ML.reactions.get_by_id("EX_cpd10515[e]").lower_bound = -0.000001
        y[7] = 0
        
    ### CYANOBACTERIA CALCULATIONS
    
    # Remove author's biomass constraint forcing which growth rate at or above experimental observation.
    # If this biomass constraint is used, all solutions at low light are infeasible.
    ML.reactions.biomass_eq_33047__hc.lower_bound = 0
    
    # Rubisco-limited photosynthesis parameters for cyanobacteria.
    ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound = -Extras.RLPC(CO2sat*4,O2sat,T)

    # Turn off HCO3 intake for cyanobacteria, just like diatom.
    ML.reactions.get_by_id("EX_cpd00242[e]").lower_bound = 0
    
    # Nitrogen fixation constraints, including boost by being symbiotic with diatom...
    ML.reactions.get_by_id("EX_cpd00528[e]").lower_bound = -(0.9*N2sat)/(165+N2sat)
    
    # If it is night, diatoms are in carbon storage mode, and photosynthesis does not occur.
    # Constraints are set accordingly.
    if Inputs['NightCheck']: 
        # Stored carbon, in the form of chrysolaminarin is consumed at night.
        # See supplemental files for our derivation of a kinetic model for chrysolaminarin consumption. 
        # Original kinetic data is from: https://doi.org/10.1016/j.ijbiomac.2023.126361
        ML.reactions.EX_chryso_e.lower_bound = -(0.935*y[4])/(133.24+y[4])
            
        # No photosynthesis occuring
        ML.reactions.EX_co2_e.lower_bound = 0
        ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound = 0
        
        # Stored carbon, in the form of chrysolaminarin is consumed at night.
        # See supplemental files for our derivation of a kinetic model for chrysolaminarin consumption. 
        # Original kinetic data is from: https://doi.org/10.1016/j.ijbiomac.2023.126361
        DML.reactions.EX_chryso_e.lower_bound = -(0.935*y[3])/(133.24+y[3]) 
        DML.reactions.EX_co2_e.lower_bound = 0
        
        # Allow a very small amount of light uptake. Not enough for photosynthesis to occur, but
        # enough to allow creation of necessary micronutrients for stored carbon metabolism.
        for k in Iphoto.keys():
            ML.reactions.get_by_id(k).upper_bound = 0
            ML.reactions.get_by_id(k).lower_bound = -1
            ML.reactions.get_by_id(k).upper_bound = -0.9999
            DML.reactions.get_by_id(k).upper_bound = 0
            DML.reactions.get_by_id(k).lower_bound = -1
            DML.reactions.get_by_id(k).upper_bound = -0.9999 
            
        # As long as silicon is available, diatoms can still grow at night using stored carbon.
        if y[5] > 0:
            ML.reactions.EX_no3_e.lower_bound = -((1.474*y[5])/(6.14+y[5]))
            ML.reactions.get_by_id("EX_cpd00209[e]").lower_bound = -((1.474*y[5])/(6.14+y[5]))
            DML.reactions.EX_no3_e.lower_bound = -((1.474*y[5])/(6.14+y[5]))
            ML.reactions.get_by_id("EX_cpd00528[e]").lower_bound = 0
    
# If it is daytime, we assume stored carbon is not consumed at all, and fixed carbon
    # is the sole food source.
    else:
        ML.reactions.EX_chryso_e.lower_bound = 0
        DML.reactions.EX_chryso_e.lower_bound = 0
        for k in Iphoto.keys():
            ML.reactions.get_by_id(k).upper_bound = 0
            ML.reactions.get_by_id(k).lower_bound = Iphoto[k]*-1.
            ML.reactions.get_by_id(k).upper_bound = Iphoto[k]*-0.9999
            DML.reactions.get_by_id(k).upper_bound = 0
            DML.reactions.get_by_id(k).lower_bound = Iphoto[k]*-1.
            DML.reactions.get_by_id(k).upper_bound = Iphoto[k]*-0.9999
        ML.reactions.EX_no3_e.lower_bound = 0
        ML.reactions.get_by_id("EX_cpd00209[e]").lower_bound = 0
    
    # The first time step requires an initial solution to be generated and saved.
    if Inputs['StartBool']:
        Inputs['StartBool'] = False
        # The simulation will always start during the daytime, startup flags accordingly.
        Inputs['NightCheck'] = False
        Inputs['FirstDay'] = True
        [Dfluxes, Inputs, DML] = fluxCalc.fluxCalc(Inputs,DML)
        [fluxes, Inputs, ML] = MergefluxCalc.MergefluxCalc(Inputs,y,T,ML,Inputs['IYSi'],Inputs['INSi'])

        refmodel.reactions.get_by_id("EX_cpd00011[e]").lower_bound = -Extras.RLPC(CO2sat*4,O2sat,T)
        Inputs['refsolution'] = refmodel.optimize()
        Inputs['LightLimit'] = Inputs['refsolution'].fluxes["EX_PHO1"]+Inputs['refsolution'].fluxes["EX_PHO2"]
        Inputs['NH4Constraint'] = -4.5*Inputs['refsolution'].fluxes["EX_cpd00528[e]"]

    # LP Calculation Flags, abort loop here if no flags are triggered and use previous LP results.
    NO3flag = (y[5] > 1e-2)
    SiO4H4flag = (y[6] > 1e-1)
    Feflag = (y[7] > 5e-5)
    
    # The first time night occurs, a basis change is necessary.
    # Calculate fluxes and flip day/night flags.
    if Inputs['NightCheck'] and Inputs['FirstNight']:
        # If enough nitrate, run the diatom-only LP, otherwise, automatically return 0 for all fluxes.
        if NO3flag:
            print("Diatom Basis Change at ", t)
            [Dfluxes, Inputs, DML] = fluxCalc.fluxCalc(Inputs,DML)
        else:
            Dfluxes['DM_biomass_c'] = 0
            Dfluxes['EX_chryso_e'] = 0
            Dfluxes['EX_o2_e'] = 0
            Dfluxes['EX_no3_e'] = 0
            Dfluxes['EX_sio4h4_e'] = 0
        # If there are enough of all limiting nutrients, run the symbiotic LP.
        # Otherwise, automatically return 0 for all fluxes.
        if NO3flag and SiO4H4flag and Feflag:
            print("Symbiosis Basis Change at ", t)
            [fluxes, Inputs, ML] = MergefluxCalc.MergefluxCalc(Inputs,y,T,ML,Inputs['IYSi'],Inputs['INSi'])
        else:
            fluxes['DM_biomass_c'] = 0
            fluxes['biomass_eq_33047__vc'] = 0
            fluxes['EX_chryso_e'] = 0
            fluxes['EX_o2_e'] = 0
            fluxes['EX_cpd00007[e]'] = 0
            fluxes['EX_no3_e'] = 0
            fluxes['EX_sio4h4_e'] = 0
            fluxes['EX_cpd00209[e]'] = 0
        Inputs['FirstNight'] = False
        Inputs['FirstDay'] = True            
    # Any other time night occurs, recalculation is likely not necessary!
    # Any necessary recalculation will be caught inside of the functions below.
    elif Inputs['NightCheck'] and not Inputs['FirstNight']:
        # If enough nitrate, run the diatom-only solver, otherwise, automatically return 0 for all fluxes.
        if NO3flag:
            [Dfluxes, Inputs, DML] = NightSolve.NightSolve(t,y,Inputs,Dfluxes,DML)
        else:
            Dfluxes['DM_biomass_c'] = 0
            Dfluxes['EX_co2_e'] = 0
            Dfluxes['EX_o2_e'] = 0
            Dfluxes['EX_no3_e'] = 0
            Dfluxes['EX_sio4h4_e'] = 0
        # If there are enough of all limiting nutrients, run the symbiotic solver.
        # Otherwise, automatically return 0 for all fluxes.    
        if NO3flag and SiO4H4flag and Feflag:
            [fluxes, Inputs, ML] = MergeNightSolve.MergeNightSolve(t,y,T,Inputs,fluxes,ML,Inputs['IYSi'],Inputs['INSi'])
        else:
            fluxes['DM_biomass_c'] = 0
            fluxes['biomass_eq_33047__vc'] = 0
            fluxes['EX_chryso_e'] = 0
            fluxes['EX_o2_e'] = 0
            fluxes['EX_cpd00007[e]'] = 0
            fluxes['EX_no3_e'] = 0
            fluxes['EX_sio4h4_e'] = 0
            fluxes['EX_cpd00209[e]'] = 0   
    # The first time day occurs, a basis change is necessary.
    # Calculate fluxes and flip day/night flags.
    elif not Inputs['NightCheck'] and Inputs['FirstDay']:
        # If enough nitrate, run the diatom-only LP, otherwise, automatically return 0 for all fluxes.
        if NO3flag:
            print("Diatom Basis Change at ", t)
            [Dfluxes, Inputs, DML] = fluxCalc.fluxCalc(Inputs,DML)
        else:
            Dfluxes['DM_biomass_c'] = 0
            Dfluxes['EX_co2_e'] = 0
            Dfluxes['EX_o2_e'] = 0
            Dfluxes['EX_no3_e'] = 0
            Dfluxes['EX_sio4h4_e'] = 0
        # If there are enough of all limiting nutrients, run the symbiotic LP.
        # Otherwise, automatically return 0 for all fluxes.
        if SiO4H4flag and Feflag:
            print("Symbiosis Basis Change at ", t)
            [fluxes, Inputs, ML] = MergefluxCalc.MergefluxCalc(Inputs,y,T,ML,Inputs['IYSi'],Inputs['INSi'])
        else:
            fluxes['DM_biomass_c'] = 0
            fluxes['biomass_eq_33047__vc'] = 0
            fluxes['EX_cpd00011[e]'] = 0
            fluxes['EX_co2_e'] = 0
            fluxes['EX_o2_e'] = 0
            fluxes['EX_cpd00007[e]'] = 0
            fluxes['EX_sio4h4_e'] = 0
            fluxes['EX_cpd10515[e]'] = 0
            fluxes['EX_cpd00528[e]'] = 0
            fluxes['carbonT'] = 0
        Inputs['FirstDay'] = False
        Inputs['FirstNight'] = True   
    # Any other time day occurs, recalculation is likely not necessary!    
    elif not Inputs['NightCheck'] and not Inputs['FirstDay']:
        # If enough nitrate, run the diatom-only solver, otherwise, automatically return 0 for all fluxes.
        if NO3flag:
            [Dfluxes, Inputs, DML] = DaySolve.DaySolve(t,y,Inputs,Dfluxes,DML)
        else:
            Dfluxes['DM_biomass_c'] = 0
            Dfluxes['EX_co2_e'] = 0
            Dfluxes['EX_o2_e'] = 0
            Dfluxes['EX_no3_e'] = 0
            Dfluxes['EX_sio4h4_e'] = 0
        # If there are enough of all limiting nutrients, run the symbiotic solver.
        # Otherwise, automatically return 0 for all fluxes.
        if SiO4H4flag and Feflag:
            [fluxes, Inputs, ML] = MergeDaySolve.MergeDaySolve(t,y,T,Inputs,fluxes,ML,Inputs['IYSi'],Inputs['INSi'])
        else:
            fluxes['DM_biomass_c'] = 0
            fluxes['biomass_eq_33047__vc'] = 0
            fluxes['EX_cpd00011[e]'] = 0
            fluxes['EX_co2_e'] = 0
            fluxes['EX_o2_e'] = 0
            fluxes['EX_cpd00007[e]'] = 0
            fluxes['EX_sio4h4_e'] = 0
            fluxes['EX_cpd10515[e]'] = 0
            fluxes['EX_cpd00528[e]'] = 0
            fluxes['carbonT'] = 0

    # In very rare cases, certain nutrients can return division by zero.
    # This usually happens at the moment when nutrient concentrations run out,
    # so zeros are returned for fluxes if this occurs.
    if math.isnan(fluxes['EX_sio4h4_e']) or math.isnan(fluxes["DM_biomass_c"]):
        return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
    # Check for nighttime! If it is night, stored carbon is being consumed by diatoms and no photosynthesis
    # is occurring. Solved fluxes are used in differential equations as below. See top of this function for
    # the order of nutrients.
    if Inputs['NightCheck']:
            
        dydt = [(Dfluxes["DM_biomass_c"]-deathK)*y[0]] # T.pseudo Biomass 0
        dydt.append((fluxes["DM_biomass_c"]-deathK)*y[1]) # C.socialis Biomass 1 
        dydt.append((fluxes['biomass_eq_33047__vc']-deathK)*y[2]) # Cyano Biomass 2
        dydt.append(Dfluxes["EX_chryso_e"]*y[0]) # T.pseudo Chyrsolaminarin 3
        dydt.append(fluxes["EX_chryso_e"]*y[1]) # C.socialis Chyrsolaminarin 4
        dydt.append(Dfluxes["EX_no3_e"]*y[0]+fluxes["EX_no3_e"]*y[1]+fluxes['EX_cpd00209[e]']*y[2]) # Nitrate 5
        dydt.append(Dfluxes["EX_sio4h4_e"]*y[0]+fluxes["EX_sio4h4_e"]*y[1]) # Silicon 6
        dydt.append(0) # Iron 7
        #dydt.append(fluxes['EX_cpd10515[e]']*y[2]) # Iron 7
        dydt.append(Dfluxes["DM_biomass_c"]*y[0]+fluxes["DM_biomass_c"]*y[1]+fluxes['biomass_eq_33047__vc']*y[2]) # Biomass without death constant 8
        dydt.append(0) # CO2 9
        dydt.append(Dfluxes["EX_o2_e"]*y[0]+fluxes["EX_o2_e"]*y[1]+fluxes['EX_cpd00007[e]']*y[2]) # O2 10
        dydt.append(0) # N2 11
            
    # Otherwise, it is daytime. If it is daytime, photosynthesis is occuring while carbon is being stored.
    # Solved fluxes are used in differential equations as below. See top of this function for
    # the order of nutrients.
    else:
        dydt = [(Dfluxes["DM_biomass_c"]-deathK)*y[0]] # T.pseudo Biomass 0
        dydt.append((fluxes["DM_biomass_c"]-deathK)*y[1]) # C.socialis Biomass 1
        dydt.append((fluxes['biomass_eq_33047__vc']-deathK)*y[2]) # Cyano Biomass
        if (Dfluxes["DM_biomass_c"]-deathK)*y[0] > 0:
            dydt.append((0.0085*math.exp(0.0163*Itot))*y[0])   # T.pseudo Chrysolaminarin
        else:
            dydt.append(0)
        if (fluxes["DM_biomass_c"]-deathK)*y[1] > 0:
            dydt.append((0.0085*math.exp(0.0163*Itot))*y[1])   # C.socialis Chyrsolaminarin 4
        else:
            dydt.append(0)
        dydt.append(Dfluxes["EX_no3_e"]*y[0]) # The symbiotic pair are forced to fix nitrogen only..
        dydt.append(Dfluxes["EX_sio4h4_e"]*y[0]+fluxes["EX_sio4h4_e"]*y[1]) # Silicon 6
        dydt.append(fluxes['EX_cpd10515[e]']*y[2]) # Iron 7
        dydt.append(Dfluxes["DM_biomass_c"]*y[0]+fluxes["DM_biomass_c"]*y[1]+fluxes['biomass_eq_33047__vc']*y[1]) # Biomass without death constant
        dydt.append(Dfluxes["EX_co2_e"]*y[0]+fluxes["EX_co2_e"]*y[1]+fluxes['EX_cpd00011[e]']*y[2]) # CO2
        dydt.append(Dfluxes["EX_o2_e"]*y[0]+fluxes["EX_o2_e"]*y[1]+fluxes['EX_cpd00007[e]']*y[2]) # O2
        dydt.append(fluxes['EX_cpd00528[e]']*y[2]) # N2
    
    return dydt