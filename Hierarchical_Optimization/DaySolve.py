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

  Outputs: 
  
        fluxes : calculated values for all fluxes relevant to system of ODEs
        Inputs : updated dictionary with new flux values and limiting reactant
        ML : same COBRA model, with possibly updated objective or constraints

  This function is very important, as it is able to detect when optimization is and is not necessary.
  For the vast majority of the simulation, CO2 will be the limiting nutrient during the day, so
  optimization rarely needs to be performed.
  
"""

def DaySolve(t,y,Inputs,Dfluxes,DML):
    
    import fluxCalc
    # Save the previous flux values in a local variable.
    prevfluxes = Dfluxes
    
    # Flag the currently limiting variable with a much higher value. 
    # Save the original value to reset later.
    saveset = DML.reactions.get_by_id(Inputs['Dprevlimit']).lower_bound
    DML.reactions.get_by_id(Inputs['Dprevlimit']).lower_bound = -100
    
    # Save array as in function fluxCalc. CO2 is used during the day instead of chryso.
    dayFluxes = ['DM_biomass_c','EX_co2_e','EX_o2_e','EX_no3_e','EX_sio4h4_e']
    
    # Save the bounds for all possibly limiting nutrients as a local variable.
    CO2bound = DML.reactions.get_by_id("EX_co2_e").lower_bound
    NO3bound = DML.reactions.get_by_id("EX_no3_e").lower_bound
    SIO4H4bound = DML.reactions.get_by_id("EX_sio4h4_e").lower_bound    
    
    # Reset the bound of the limiting reactant to its original value.
    DML.reactions.get_by_id(Inputs['Dprevlimit']).lower_bound = saveset
    fluxval = list(Dfluxes.values())
    
    # These determine if the saved solution's limiting nutrient is still the limiting reactant.
    # If any of these come back false, the basis set needs to be updated via optimization.
    CO2bool = (CO2bound < fluxval[1])
    NO3bool = (NO3bound < fluxval[3])
    SIO4H4bool = (SIO4H4bound < fluxval[4])
    
    if CO2bool and NO3bool and SIO4H4bool:
        Dfluxes = {}
        # The system is linear. If the limiting nutrient hasn't changed, one can simply
        # scale all fluxes by the change in the limiting nutrient's constraint.
        for idx, name in enumerate(prevfluxes):
            Dfluxes[dayFluxes[idx]] = Inputs['Dprevflux'][name]*(saveset/Inputs['Dsavedflux'])
    else:
        # Perform the optimization and save new fluxes.
        print("Diatom Basis Change at ", t)
        [Dfluxes, Inputs, DML] = fluxCalc.fluxCalc(Inputs, DML)
        Inputs['Dprevflux'] = Dfluxes
        
    return Dfluxes, Inputs, DML