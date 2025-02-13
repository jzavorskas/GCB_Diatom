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
        DML : COBRA model structure for Thalassiosira pseudonana containing all current constraints

  Outputs: 
  
        fluxes : calculated values for all fluxes relevant to system of ODEs
        Inputs : updated dictionary with new flux values and limiting reactant
        DML : same COBRA model, with possibly updated objective or constraints

  This function is very important, as it is able to detect when optimization is and is not necessary.
  For the vast majority of the simulation, chrysolaminarin will be the limiting nutrient at night, so
  optimization rarely needs to be performed.
  
"""

def NightSolve(t,y,Inputs,Dfluxes,DML):
    
    import fluxCalc
    # Save the previous flux values in a local variable.
    prevfluxes = Dfluxes
    
    # Flag the currently limiting variable with a much higher value. 
    # Save the original value to reset later.
    saveset = DML.reactions.get_by_id(Inputs['Dprevlimit']).lower_bound
    DML.reactions.get_by_id(Inputs['Dprevlimit']).lower_bound = -100
    
    # Save array as in function fluxCalc. Chryso is used at night instead of CO2.
    nightFluxes = ['DM_biomass_c', 'EX_chryso_e','EX_o2_e','EX_no3_e','EX_sio4h4_e']
    
    # Save the bounds for all possibly limiting nutrients as a local variable.
    Chrysobound = DML.reactions.get_by_id("EX_chryso_e").lower_bound
    NO3bound = DML.reactions.get_by_id("EX_no3_e").lower_bound
    SIO4H4bound = DML.reactions.get_by_id("EX_sio4h4_e").lower_bound    
    
    # Reset the bound of the limiting reactant to its original value.
    DML.reactions.get_by_id(Inputs['Dprevlimit']).lower_bound = saveset
    fluxval = list(Dfluxes.values())
    
    # These determine if the saved solution's limiting nutrient is still the limiting reactant.
    # If any of these come back false, the basis set needs to be updated via optimization.
    Chrysobool = (Chrysobound - fluxval[1]) < -0.001
    NO3bool = (NO3bound - fluxval[3]) < -0.001
    SIO4H4bool = (SIO4H4bound - fluxval[4]) < -0.001

    if Chrysobool and NO3bool and SIO4H4bool:
        Dfluxes = {}
        # The system is linear. If the limiting nutrient hasn't changed, one can simply
        # scale all fluxes by the change in the limiting nutrient's constraint.
        for idx, name in enumerate(prevfluxes):
            Dfluxes[nightFluxes[idx]] = Inputs['Dprevflux'][name]*(saveset/Inputs['Dsavedflux'])
    else:
        # Perform the optimization and save new fluxes.
        print("Diatom Basis Change at ", t)
        [Dfluxes, Inputs, DML] = fluxCalc.fluxCalc(Inputs, DML)
        Inputs['Dprevflux'] = Dfluxes
        
    return Dfluxes, Inputs, DML