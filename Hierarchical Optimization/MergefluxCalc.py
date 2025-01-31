"""
_____________________________________________________________
Function: Flux Calculation (fluxCalc) function

  Inputs: 
  
        Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations
        ML : COBRA model structure for Thalassiosira pseudonana containing all current constraints

  Outputs: 
  
        fluxes : calculated values for all fluxes relevant to system of ODEs
        Inputs : updated dictionary with new flux values and limiting reactant
        ML : same COBRA model, with possibly updated objective or constraints

  As addressed in the text, the optimal cost of the LP is unique but the other fluxes in the solution 
  may not be. As such, to guarantee unique solutions for all fluxes used to update the ODEs, one must
  perform a hierarchical optimization. Specifically, all fluxes used in the ODE should be made the objective
  function one-by-one, while adding the previous optimal fluxes as constraints on the problem. In the case 
  of this function, once biomass is calculated, the optimal value would be fixed for all subsequent
  optimizations.
  
"""

def MergefluxCalc(Inputs,y,T,ML,IYSi,INSi):
    
    import MergeTransfer, Extras
    import cobra
    
    ML.reactions.EX_co2_e.lower_bound = -Extras.RLPD(Inputs['sat'][0]*4,Inputs['sat'][1],T)
    
    # Set biomass maximization as objective.
    ML.objective = {ML.reactions.get_by_id("DM_biomass_c"):1,ML.reactions.get_by_id("biomass_eq_33047__vc"):1}
    # Initialize dictionary to hold fluxes with their reaction IDs as keys.
    fluxes = {}
    
    # Initialize list to hold constraints.
    Cons = []
    
    # Check if it is daytime...
    if not Inputs['NightCheck']:
        
        # During the day, CO2 is used for carbon fixation, so it is the second in hierarchical optimization.
        inFluxes = ['DM_biomass_c','biomass_eq_33047__vc','EX_cpd00011[e]','EX_co2_e',
                    'EX_o2_e','EX_cpd00007[e]',
                    'EX_sio4h4_e','EX_cpd10515[e]','EX_cpd00528[e]','carbonT']
        dyFluxes = ['EX_cpd00011[e]','EX_co2_e',
                    'EX_o2_e','EX_cpd00007[e]',
                    'EX_sio4h4_e','EX_cpd10515[e]','EX_cpd00528[e]','carbonT']
        
        # Other than biomass, all other optimizations are performed as minimizations to avoid unnecessary
        # uptake for peripheral pathways that wouldn't happen during normal growth (i.e. don't show up if 
        # biomass is maximized).
        optsign = [1,1,-1,-1,-1,1,1,-1] #CyanoCO2,DiatomCO2,DiatomO2,CyanoO2,DiatomSi(OH)4,CyanoFe2+,CyanoN2,CarbonTransfer
        
        # Add constraint for light concentration...
        Cons.append(ML.problem.Constraint(
                (ML.reactions.EX_PHO1.flux_expression + ML.reactions.EX_PHO2.flux_expression),
                lb=-Inputs['Itot'],
                ub=0))
        ML.add_cons_vars(Cons[0])

        # Add constraint for nitrogen transfer...
        Cons.append(ML.problem.Constraint(
                ML.reactions.nh4T.flux_expression,
                lb=0,
                ub=Inputs['NH4Constraint']))
        ML.add_cons_vars(Cons[1])

        # Add constraint for carbon transfer...
        if not Inputs['NightCheck']:
            [CarbonConstraint, CT] = MergeTransfer.MergeTransfer(Inputs,y,T,ML,IYSi,INSi)
            Cons.append(CarbonConstraint)
        # Unless it is nighttime, then just force biomass expressions to be equal.
        else:
            Cons.append(ML.problem.Constraint(
                (ML.reactions.DM_biomass_c.flux_expression - ML.reactions.biomass_eq_33047__vc.flux_expression),
                lb=0,
                ub=0))   

        ML.add_cons_vars(Cons[2])  
    
    # Or if it is night! 
    else:
        ML.reactions.EX_co2_e.lower_bound = 0
        ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound = 0
        
        # At night, chrysolaminarin is consumed autotrophically, so it is second.
        inFluxes = ['DM_biomass_c','biomass_eq_33047__vc', 'EX_chryso_e',
                    'EX_o2_e','EX_cpd00007[e]', 'EX_no3_e','EX_sio4h4_e',
                    'EX_cpd00209[e]']
        dyFluxes = ['EX_chryso_e','EX_o2_e','EX_cpd00007[e]', 'EX_no3_e',
                    'EX_sio4h4_e','EX_cpd00209[e]']
        
        # Other than biomass, all other optimizations are performed as minimizations to avoid unnecessary
        # uptake for peripheral pathways that wouldn't happen during normal growth (i.e. don't show up if 
        # biomass is maximized).
        optsign = [1,-1,1,1,-1,1,-1] #DiatomChryso,DiatomO2,CyanoO2,DiatomNO3,CyanoNO3,DiatomSi(OH)4,CyanoFe2+,CarbonTransfer     
        
        Cons.append(ML.problem.Constraint(
            (ML.reactions.EX_PHO1.flux_expression + ML.reactions.EX_PHO2.flux_expression),
            lb=-1,
            ub=0))
        ML.add_cons_vars(Cons[0])

        Cons.append(ML.problem.Constraint(
            (ML.reactions.DM_biomass_c.flux_expression - ML.reactions.biomass_eq_33047__vc.flux_expression),
            lb=0,
            ub=0))
        ML.add_cons_vars(Cons[1])

        Cons.append(ML.problem.Constraint(
            (ML.reactions.DM_biomass_c.flux_expression + ML.reactions.biomass_eq_33047__vc.flux_expression),
            lb=5e-7,
            ub=1000))
        ML.add_cons_vars(Cons[2])
    
    # Optimize to maximize biomass, save the value.
    solution = ML.optimize()
    fluxes['DM_biomass_c'] = solution.fluxes['DM_biomass_c']
    fluxes['biomass_eq_33047__vc'] = solution.fluxes['biomass_eq_33047__vc']

    Cons.append(ML.problem.Constraint(
        getattr(ML.reactions, 'DM_biomass_c').flux_expression,
        lb=fluxes['DM_biomass_c'],
        ub=fluxes['DM_biomass_c']))
    ML.add_cons_vars(Cons[3])
                
    Cons.append(ML.problem.Constraint(
        getattr(ML.reactions, 'biomass_eq_33047__vc').flux_expression,
        lb=fluxes['biomass_eq_33047__vc'],
        ub=fluxes['biomass_eq_33047__vc']))
    ML.add_cons_vars(Cons[4])
                
    for idx, dyFlux in enumerate(dyFluxes):
        
        # Set the current flux in the loop as the objective function, its optimization sign
        # will make the objective a minimization.
        ML.objective = {getattr(ML.reactions, dyFlux):optsign[idx]}

        # Optimize the current flux and save its value in the flux dictionary.
        solution = ML.optimize()
        fluxes[dyFlux] = solution.fluxes[dyFlux]        
              
        # Set the previously calculated flux as a new constraint.
        Cons.append(ML.problem.Constraint(
            getattr(ML.reactions, dyFlux).flux_expression,
            lb=fluxes[dyFlux],
            ub=fluxes[dyFlux]))
        ML.add_cons_vars(Cons[idx+5])

    # VERY IMPORTANT: Remove all constraints before moving on!
    for removeCons in Cons:
        ML.remove_cons_vars(removeCons)
    
    # Find the limiting flux (which flux has the same value as its constraint?)
    # Save the limiting flux's name and value as separate entries in Inputs.
    for idx, flux in enumerate(fluxes.values()):
        # Don't check biomass values...
        if idx <= 1:
            pass
        elif (abs(flux - (getattr(ML.reactions, inFluxes[idx]).lower_bound))) <= 0.000001: 
            Inputs['prevlimit'] = inFluxes[idx]
            Inputs['savedflux'] = flux
    
    # Save the current flux profile in inputs, to use later.        
    Inputs['prevflux'] = fluxes
    return fluxes, Inputs, ML