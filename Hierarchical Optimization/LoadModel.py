"""
_____________________________________________________________
Function: Diatom Model Loading (LoadModel) function

  Inputs: 
  
        NONE

  Outputs: 
  
        ML : COBRA model for Thalassiosira pseudonana, with default objective and constraints

  As addressed in the text, the optimal cost of the LP is unique but the other fluxes in the solution 
  may not be. As such, to guarantee unique solutions for all fluxes used to update the ODEs, one must
  perform a hierarchical optimization. Specifically, all fluxes used in the ODE should be made the objective
  function one-by-one, while adding the previous optimal fluxes as constraints on the problem. In the case 
  of this function, once biomass is calculated, the optimal value would be fixed for all subsequent
  optimizations.
  
"""

def LoadModel(Inputs): 
    import cobra
    # Import iTps1432 Thalassiosira pseudonana
    ML = cobra.io.read_sbml_model('Models/Thaps_ML_model_n.xml')
    ML.objective = {ML.reactions.get_by_id("DM_biomass_c"):1}

    # Create reaction objects for each of the new reactions we are importing.
    # The first reaction is created implicitly below as an exchange flux for stored chryso.
    # The two here are a "transport" reaction, unstoring the chryso, and a conversion directly to glucose.
    gluconv = cobra.Reaction('gluconv')
    chrysoT = cobra.Reaction('chrysoT')
    GLCt1 = cobra.Reaction('GLCt1')
    G6PD_c = cobra.Reaction('G6PD_c')
    ML.add_reactions([G6PD_c])
    ML.reactions.G6PD_c.build_reaction_from_string('g6p_c --> glc__D_c')
    
    # Create metabolites to represent stored (represented as "extracellular") and available chryso.
    # Also add missing glucose extracellular metabolite.
    ML.add_metabolites([
        cobra.Metabolite(
        'chryso_c',
        name='chrysolaminarin',
        compartment='c'
        ),
        cobra.Metabolite(
        'chryso_e',
        name='chrysolaminarin',
        compartment='e'
        ),
        cobra.Metabolite(
        'glc__D_e',
        name='D-Glucose',
        compartment='e'
        ),
    ])

    # Add new reaction items and associated exchange reactions,
    # fill added reactions with reaction expressions.
    ML.add_reactions([GLCt1,chrysoT,gluconv])
    ML.add_boundary(ML.metabolites.get_by_id("glc__D_e"), type="exchange")
    ML.add_boundary(ML.metabolites.get_by_id("chryso_e"), type="exchange")
    ML.reactions.gluconv.build_reaction_from_string('chryso_c --> glc__D_c')
    ML.reactions.chrysoT.build_reaction_from_string('chryso_e --> chryso_c')
    ML.reactions.GLCt1.build_reaction_from_string('glc__D_e --> glc__D_c')

    # Must run this cell or it will rely entirely on heterotrophic growth!!!
    # Defaults for co2 exchange are ub = 0, lb = 0.
    ML.reactions.EX_co2_e.upper_bound = 1000
    ML.reactions.EX_co2_e.lower_bound = -1000
    
    # Run this cell if you wish to perform autotrophic growth only.
    ML.reactions.EX_glc__D_e.lower_bound = 0
    ML.reactions.EX_glc__D_e.upper_bound = 0
    
    # Shut down or turn on the chrysolaminarin reaction.
    ML.reactions.EX_chryso_e.lower_bound = 0
    
    # Default constraint values, which will likely get changed in the dFBA loop.
    # If not involve in the ODE system, constraints are left open.
    ML.reactions.EX_no3_e.lower_bound = -0.1
    ML.reactions.EX_no2_e.lower_bound = 0
    ML.reactions.EX_nh4_e.lower_bound = 0
    ML.reactions.EX_pi_e.lower_bound = -1000
    ML.reactions.EX_so4_e.lower_bound = -1
    ML.reactions.EX_sio4h4_e.lower_bound = -1000
    ML.reactions.EX_cncbl3_e.lower_bound = -1000
    ML.reactions.EX_co2_e.lower_bound = -10
    ML.reactions.EX_hco3_e.lower_bound = 0
    
    return ML