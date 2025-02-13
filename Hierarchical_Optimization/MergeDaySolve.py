"""
_____________________________________________________________
Function: Day Solving (DaySolve) function

  Inputs: 
  
        t : Current ODE solver time
        y : Current ODE solver y values (biomass, nutrients, etc.)
        Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations
        fluxes : current fluxes, either calculated via hierarchical optimization or held over from a 
                 previous one
        ML : COBRA model structure for Thalassiosira pseudonana containing all current constraints
        IYSi : reference table when silicon is NOT limiting
        INSi : reference table when silicon IS limiting

  Outputs: 
  
        fluxes : calculated values for all fluxes relevant to system of ODEs
        Inputs : updated dictionary with new flux values and limiting reactant
        ML : same COBRA model, with possibly updated objective or constraints

  This function is very important, as it is able to detect when optimization is and is not necessary.
  For the vast majority of the simulation, CO2 will be the limiting nutrient during the day, so
  optimization rarely needs to be performed.
  
"""

def MergeDaySolve(t,y,T,Inputs,fluxes,ML,IYSi,INSi):
    
    import MergeTransfer, MergefluxCalc
    
    # Save the previous flux values in a local variable.
    prevfluxes = fluxes 
    
    # Flag the currently limiting variable with a much higher value. 
    # Save the original value to reset later.
    saveset = ML.reactions.get_by_id(Inputs['prevlimit']).lower_bound
    ML.reactions.get_by_id(Inputs['prevlimit']).lower_bound = -100
    
    # Save array as in function fluxCalc. CO2 is used during the day instead of chryso.
    dayFluxes = ['DM_biomass_c','biomass_eq_33047__vc','EX_cpd00011[e]','EX_co2_e',
                    'EX_o2_e','EX_cpd00007[e]',
                    'EX_sio4h4_e','EX_cpd10515[e]','EX_cpd00528[e]','carbonT']
    
    # Save the bounds for all possibly limiting nutrients as a local variable.
    DCO2bound = ML.reactions.get_by_id("EX_co2_e").lower_bound
    CCO2bound = ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound
    SIO4H4bound = ML.reactions.get_by_id("EX_sio4h4_e").lower_bound
    CFebound = ML.reactions.get_by_id("EX_cpd10515[e]").lower_bound
    CN2bound = ML.reactions.get_by_id("EX_cpd00528[e]").lower_bound
    
    # Reset the bound of the limiting reactant to its original value.
    ML.reactions.get_by_id(Inputs['prevlimit']).lower_bound = saveset
    fluxval = list(fluxes.values())
    
    DCO2boundreal = ML.reactions.get_by_id("EX_co2_e").lower_bound
    CCO2boundreal = ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound
    
    # These determine if the saved solution's limiting nutrient is still the limiting reactant.
    # If any of these come back false, the basis set needs to be updated via optimization.
    DCO2bool = (DCO2bound <= fluxval[3])
    CCO2bool = (CCO2bound <= fluxval[2])
    SIO4H4bool = (SIO4H4bound <= fluxval[6])
    CN2bool = (CN2bound <= fluxval[8])
    
    Megabool = (DCO2bool and CCO2bool and SIO4H4bool and CN2bool)
    
    if Megabool:
        
        if Inputs['prevlimit'] == "EX_co2_e":
            
            DCO2prev = prevfluxes['EX_co2_e']
            CCO2prev = prevfluxes['EX_cpd00011[e]']
            [CarbonTransfer, CT] = MergeTransfer.MergeTransfer(Inputs,y,T,ML,IYSi,INSi)
            
            CO2boundsum = DCO2boundreal + CCO2boundreal
            CO2prevsum = DCO2prev + CCO2prev
                 
            fluxes = {}
            
            fluxes['DM_biomass_c'] = prevfluxes['DM_biomass_c']*(CO2boundsum/CO2prevsum)
            fluxes['biomass_eq_33047__vc'] = prevfluxes['biomass_eq_33047__vc']*(CO2boundsum/CO2prevsum)
            fluxes['EX_cpd00011[e]'] = prevfluxes['EX_cpd00011[e]']*(CCO2boundreal/CCO2prev)
            fluxes['EX_co2_e'] = prevfluxes['EX_co2_e']*(DCO2boundreal/DCO2prev)
            fluxes['EX_o2_e'] = prevfluxes['EX_o2_e']*(DCO2boundreal/DCO2prev)
            fluxes['EX_cpd00007[e]'] = prevfluxes['EX_cpd00007[e]']*((CT*1.2855+CCO2bound)/(prevfluxes['carbonT']*1.2855+prevfluxes['EX_cpd00011[e]']))
            fluxes['EX_sio4h4_e'] = prevfluxes['EX_sio4h4_e']*(CO2boundsum/CO2prevsum)
            fluxes['EX_cpd10515[e]'] = prevfluxes['EX_cpd10515[e]']*(CO2boundsum/CO2prevsum)
            fluxes['EX_cpd00528[e]'] = prevfluxes['EX_cpd00528[e]']*(CO2boundsum/CO2prevsum)
            fluxes['carbonT'] = CT      
                
        else:
            
            print("Non-CO2 Limiter - Running optimization")
            [fluxes, Inputs, ML] = MergefluxCalc.MergefluxCalc(Inputs,y,T,ML,IYSi,INSi)
            Inputs['prevflux'] = fluxes
            
    else:
        # Perform the optimization and save new fluxes.
        print("Symbiosis Basis Change at ", t)
        [fluxes, Inputs, ML] = MergefluxCalc.MergefluxCalc(Inputs,y,T,ML,IYSi,INSi)
        Inputs['prevflux'] = fluxes
        
    return fluxes, Inputs, ML