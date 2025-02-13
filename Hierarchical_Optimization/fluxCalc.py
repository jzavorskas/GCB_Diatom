"""
_____________________________________________________________
Function: Flux Calculation (fluxCalc) function

  Inputs: 
  
        Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations
        DML : COBRA model structure for Thalassiosira pseudonana containing all current constraints

  Outputs: 
  
        fluxes : calculated values for all fluxes relevant to system of ODEs
        Inputs : updated dictionary with new flux values and limiting reactant
        DML : same COBRA model, with possibly updated objective or constraints

  As addressed in the text, the optimal cost of the LP is unique but the other fluxes in the solution 
  may not be. As such, to guarantee unique solutions for all fluxes used to update the ODEs, one must
  perform a hierarchical optimization. Specifically, all fluxes used in the ODE should be made the objective
  function one-by-one, while adding the previous optimal fluxes as constraints on the problem. In the case 
  of this function, once biomass is calculated, the optimal value would be fixed for all subsequent
  optimizations.
  
"""

def fluxCalc(Inputs, DML):
    
    # Set biomass maximization as objective.
    DML.objective = {getattr(DML.reactions, "DM_biomass_c"):1}
    # Initialize dictionary to hold fluxes with their reaction IDs as keys.
    Dfluxes = {}
    
    # Check if it is daytime...
    if not Inputs['NightCheck']:
        
        # During the day, CO2 is used for carbon fixation, so it is the second in hierarchical optimization.
        inFluxes = ['DM_biomass_c','EX_co2_e','EX_o2_e','EX_no3_e','EX_sio4h4_e']
        dyFluxes = ['EX_co2_e','EX_o2_e','EX_no3_e','EX_sio4h4_e']
        
        # Other than biomass, all other optimizations are performed as minimizations to avoid unnecessary
        # uptake for peripheral pathways that wouldn't happen during normal growth (i.e. don't show up if 
        # biomass is maximized).
        optsign = [1,-1,1,1]
        
    # Or if it is night! 
    else:
        # At night, chrysolaminarin is consumed autotrophically, so it is second.
        inFluxes = ['DM_biomass_c', 'EX_chryso_e','EX_o2_e','EX_no3_e','EX_sio4h4_e']
        dyFluxes = ['EX_chryso_e','EX_o2_e','EX_no3_e','EX_sio4h4_e']
        
        # Other than biomass, all other optimizations are performed as minimizations to avoid unnecessary
        # uptake for peripheral pathways that wouldn't happen during normal growth (i.e. don't show up if 
        # biomass is maximized).
        optsign = [1,-1,1,1]      
        
    # Initialize list to hold constraints.
    Cons = []

    # Optimize to maximize biomass, save the value.
    solution = DML.optimize()
    Dfluxes[inFluxes[0]] = solution.fluxes['DM_biomass_c']

    
    for idx, dyFlux in enumerate(dyFluxes):
        
        # Set the current flux in the loop as the objective function, its optimization sign
        # will make the objective a minimization.
        DML.objective = {getattr(DML.reactions, dyFlux):optsign[idx]}

        # Set the previously calculated flux as a new constraint.
        Cons.append(DML.problem.Constraint(
            getattr(DML.reactions, inFluxes[idx]).flux_expression,
            lb=Dfluxes[inFluxes[idx]],
            ub=Dfluxes[inFluxes[idx]]))
        DML.add_cons_vars(Cons[idx])

        # Optimize the current flux and save its value in the flux dictionary.
        solution = DML.optimize()
        Dfluxes[inFluxes[idx+1]] = solution.fluxes[dyFlux]

    # VERY IMPORTANT: Remove all constraints before moving on!
    for removeCons in Cons:

        DML.remove_cons_vars(removeCons)
    
    # Find the limiting flux (which flux has the same value as its constraint?)
    # Save the limiting flux's name and value as separate entries in Inputs.
    for idx, flux in enumerate(Dfluxes.values()):
        # Don't check biomass...
        if idx == 0:
            pass
        elif (abs(flux - (getattr(DML.reactions, inFluxes[idx]).lower_bound))) <= 0.0001: 
            Inputs['Dprevlimit'] = inFluxes[idx]
            Inputs['Dsavedflux'] = flux
            
    # Save the current flux profile in inputs, to use later.        
    Inputs['Dprevflux'] = Dfluxes
    return Dfluxes, Inputs, DML