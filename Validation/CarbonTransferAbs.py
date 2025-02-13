"""
Author: Joe Zavorskas
Date Started: 04/11/2024
Last Edit: 02/11/2025
____________________________________________________________________
Function: Absolute Carbon Transfer (CarbonTransferAbs) calculation

  Inputs: 
  
        T : Temperature at current loop iteration (degrees Celcius)
        NO3 : Nitrate concentration at current loop iteration (umol/L)
        Sili : Silicic acid concentration at current loop iteration (umol/L)
        ML : COBRA model structure for the symbiotic pairing:
                Thalassiosira pseudonana and Anabaena variabilis containing all current constraints
        refmodel : COBRA model structure for Anabaena variabilis ONLY containing all current constraints;
                    to be used as a reference for the cyanobacteria's basal metabolism
        Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations

  Outputs: 
  
        CarbonTransferVal : The calculated carbon transfer value, given the current loop's value for
                            temperature, nitrate concentration, and silicic acid concentration

This file contains a calculation function which is used in tandem with "CarbonTransferTable.py."
Together, these files generate a lookup table of the carbon transfer value that leads to equal 
growth rates for diatom and cyanobacteria at a range of temperatures, nitrate concentrations, and
silicic acid concentrations.

"""

def CarbonTransferAbs(T,NO3,Sili,ML,refmodel,Inputs):

    import Extras, MergeLoadModel

    Henry = Extras.henry(T)
    CO2sat = Henry[0]*(4.21e-4) # Multiply Henry coefficient by: Partial pressure of CO2 in atm.
    O2sat = Henry[1]*(.2105) # Partial pressure of O2 in atm
    N2sat = Henry[2]*(.781) # Partial pressure of N2 in atm

    # This block solves a "reference model", which will give Anabaena's basal consumption of CO2
    # and generation of fixed nitrogen from N2.
    refmodel.reactions.get_by_id("EX_cpd00011[e]").lower_bound = -Extras.RLPC(CO2sat*4,O2sat,T)
    refmodel.reactions.get_by_id("EX_cpd00242[e]").lower_bound = 0
    refmodel.reactions.get_by_id("EX_cpd00027[e]").lower_bound = 0
    refmodel.reactions.get_by_id("EX_cpd00528[e]").lower_bound = -(0.9*N2sat)/(165+N2sat)
    refmodel.reactions.biomass_eq_33047__hc.lower_bound = 0
    refsolution = refmodel.optimize()

    # This block addresses infeasibility issues when concentrations are *very low*.
    # Slightly loosens toleranaces and institutes a quick timeout if solver is stuck.
#    ML.solver.configuration.verbosity = 3
    ML.solver.configuration.timeout = 5
    ML.solver.configuration.tolerance_feasibility = 1e-5


    # Use Temp-dependent Rubisco enzyme kinetics to set bounds of CO2 consumption for diatom.
    # Data can be found at:
    #   1) https://doi.org/10.1007/s11120-014-0067-8
    #   2) https://doi.org/10.1093/jxb/erw267
    ML.reactions.EX_co2_e.upper_bound = 1000
    ML.reactions.EX_co2_e.lower_bound = -Extras.RLPD(CO2sat*4,O2sat,T)
    ML.reactions.EX_hco3_e.lower_bound = 0

    # Michaelis-Menten expression for NO3 consumption by diatom. Data can be found at:
    # https://doi.org/10.3354/meps09088
    ML.reactions.EX_no3_e.lower_bound = -((1.474*NO3)/(6.14+NO3))

    # Hill exponent for cooperative allosteric interactions in the silicon incorporation process.
    # Modified Michaelis-Menten with hill kinetics for uptake of silicic acid by diatom:
    # Data can be found at: https://doi.org/10.1104/pp.107.107094
    n = 1.9
    ML.reactions.EX_sio4h4_e.lower_bound = -((1.961*(Sili**n))/((8.1**n)+(Sili**n)))

    # Use Temp-dependent Rubisco enzyme kinetics to set bounds of CO2 consumption for cyanobacteria.
    # Data can be found at the same links as diatom above.
    ML.reactions.get_by_id("EX_cpd00011[e]").upper_bound = 1000
    ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound = -Extras.RLPC(CO2sat*4,O2sat,T)
    # Force cyanobacteria to produce CO2 at the exact maximum possible rate.
    ML.reactions.get_by_id("EX_cpd00011[e]").upper_bound = -Extras.RLPC(CO2sat*4,O2sat,T)*0.9999

    # Michaelis-Menten saturation kinetics for cyanobacteria nitrogen fixation shown below.
    # Data can be found at: https://doi.org/10.1111/j.1574-6968.1981.tb07607.x
    # Parameters were modified based on knowledge that cyanobacteria in symbiosis reach 900%
    # of their basal nitrogen fixation rate. We have conservatively multiply the maximum
    # turnover rate by 2.5 to compensate. Symbiosis data here: https://doi.org/10.3390/plants9020192
    ML.reactions.get_by_id("EX_cpd00528[e]").lower_bound = -(0.9*N2sat)/(165+N2sat)

    # Add a constraint to the FBA linear program that allows the cyanobacteria to uptake
    # only up to the amount of available light.
    Imax = ML.problem.Constraint(
        (ML.reactions.EX_PHO1.flux_expression + ML.reactions.EX_PHO2.flux_expression),
        lb=-Inputs['Itot'],
        ub=0)
    ML.add_cons_vars(Imax)
    
    # Add a constraint to the FBA linear program that requires the cyanobacterium to transfer
    # exactly 450% of its basal nitrogen fixation rate to the diatom. We take a conservative
    # estimate here again, having calculated that at 5 degreesC, cyanobacteria transfer 550% of
    # their basal nitrogen rate.
    AmmoniaConstraint = ML.problem.Constraint(
        ML.reactions.nh4T.flux_expression,
        lb=-4.5*refsolution.fluxes["EX_cpd00528[e]"],
        ub=-4.5*refsolution.fluxes["EX_cpd00528[e]"])
    ML.add_cons_vars(AmmoniaConstraint)
    
    # Add a constraint to the FBA linear program that forces the growth rates of diatom and
    # cyanobacteria to be equal. This is the most important part of this calculation.
    BioConstraint = ML.problem.Constraint(
        (ML.reactions.DM_biomass_c.flux_expression - ML.reactions.biomass_eq_33047__vc.flux_expression),
        lb=0,
        ub=0)
    ML.add_cons_vars(BioConstraint)
    
    ### In very rare situations, the FBA linear program will respond to the above code block by 
    ### making both biomasses zero. If this occurs, use the following constraint to force a
    ### non-trivial solution.
 #   NonTrivial =  ML.problem.Constraint(
  #      (ML.reactions.DM_biomass_c.flux_expression + ML.reactions.biomass_eq_33047__vc.flux_expression),
   #     lb=0.005,
    #    ub=1000)
    #ML.add_cons_vars(NonTrivial)
    
    # Solve and store the FBA problem.
    solution = ML.optimize()
    
    # It is very important to remove all constraints before carrying out the next loop,
    # or difficult to debug infeasibility will arise.
    ML.remove_cons_vars(Imax)
    ML.remove_cons_vars(AmmoniaConstraint)
    ML.remove_cons_vars(BioConstraint)
   # ML.remove_cons_vars(NonTrivial)

    # User-facing, will print out all the fluxes calculated on this loop.
    print(solution.fluxes["carbonT"])
    print(ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound)
    print(ML.reactions.EX_co2_e.lower_bound)
    print(solution.fluxes["DM_biomass_c"])
    print(solution.fluxes["biomass_eq_33047__vc"])
    print(solution.fluxes["EX_cpd00011[e]"])
    
    # Ultimately, the carbon transfer is calculated as the biological flux of individual carbon atoms
    # (originally reported as transfer of glucose in "carbonT"), divided by basal cyanobacteria
    # carbon production. This value will eventually tell the FBA program how much of the cyanobacteria's
    # metabolism needs to be supplemented by the diatom to make growth rates equal.
    CarbonTransferVal = (6*solution.fluxes["carbonT"]/(-ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound))

    return CarbonTransferVal