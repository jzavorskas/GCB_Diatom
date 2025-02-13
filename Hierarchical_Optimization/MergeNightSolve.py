"""
_____________________________________________________________
Function: Night Solving (NightSolve) function

  Inputs: 
  
        t : Current ODE solver time
        y : Current ODE solver y values (biomass, nutrients, etc.)
        Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations
        fluxes : current fluxes, either calculated via hierarchical optimization or held over from a 
                 previous one
        ML : COBRA model structure for Thalassiosira pseudonana containing all current constraints

  Outputs: 
  
        fluxes : calculated values for all fluxes relevant to system of ODEs
        Inputs : updated dictionary with new flux values and limiting reactant
        ML : same COBRA model, with possibly updated objective or constraints

  This function is very important, as it is able to detect when optimization is and is not necessary.
  For the vast majority of the simulation, chrysolaminarin will be the limiting nutrient at night, so
  optimization rarely needs to be performed.
  
"""

def MergeNightSolve(t,y,T,Inputs,fluxes,ML,IYSi,INSi):
    
    import MergeTransfer, MergefluxCalc
    
    # Save the previous flux values in a local variable.
    prevfluxes = fluxes
    
    # Flag the currently limiting variable with a much higher value. 
    # Save the original value to reset later.
    saveset = ML.reactions.get_by_id(Inputs['prevlimit']).lower_bound
    ML.reactions.get_by_id(Inputs['prevlimit']).lower_bound = -100
    
    # Save array as in function fluxCalc. Chryso is used at night instead of CO2.
    nightFluxes = ['DM_biomass_c','biomass_eq_33047__vc', 'EX_chryso_e',
                    'EX_o2_e','EX_cpd00007[e]','EX_no3_e','EX_sio4h4_e','EX_cpd00209[e]']
    
    # Save the bounds for all possibly limiting nutrients as a local variable.
    DChrysobound = ML.reactions.get_by_id("EX_chryso_e").lower_bound
    DNO3bound = ML.reactions.get_by_id("EX_no3_e").lower_bound
    CNO3bound = ML.reactions.get_by_id("EX_cpd00209[e]").lower_bound    
    DSIO4H4bound = ML.reactions.get_by_id("EX_sio4h4_e").lower_bound
    CFebound = ML.reactions.get_by_id("EX_cpd10515[e]").lower_bound   
    
    # Reset the bound of the limiting reactant to its original value.
    ML.reactions.get_by_id(Inputs['prevlimit']).lower_bound = saveset
    fluxval = list(fluxes.values())
    
    # These determine if the saved solution's limiting nutrient is still the limiting reactant.
    # If any of these come back false, the basis set needs to be updated via optimization.
    DChrysobool = (DChrysobound <= fluxval[2]) or abs(DChrysobound - fluxval[2]) <= 1e3
    DSIO4H4bool = (DSIO4H4bound <= fluxval[6]) 
    CNO3bool = (CNO3bound <= fluxval[7]) 
    
    Megabool = (DChrysobool and DSIO4H4bool and CNO3bool)
     
    if Megabool:
        fluxes = {}
        # The system is linear. If the limiting nutrient hasn't changed, one can simply
        # scale all fluxes by the change in the limiting nutrient's constraint.
        for idx, name in enumerate(nightFluxes):
            fluxes[nightFluxes[idx]] = Inputs['prevflux'][name]*(saveset/Inputs['savedflux'])
    else:
        # Perform the optimization and save new fluxes.
        print("Symbiosis Basis Change at ", t)
        [fluxes, Inputs, ML] = MergefluxCalc.MergefluxCalc(Inputs,y,T,ML,IYSi,INSi)
        Inputs['prevflux'] = fluxes
        
    fluxes["EX_co2_e"] = 0
    fluxes["EX_cpd00011[e]"] = 0
        
    return fluxes, Inputs, ML