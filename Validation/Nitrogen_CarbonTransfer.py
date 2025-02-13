"""
_________________________________________________________________________________________
Function: nitrogen and carbon transfer validation (Nitrogen_CarbonTransfer) calculation

  Inputs: 
  
        Inputs : a dictionary containing datetime information for light calculation
                 and information about previous flux calculations
        T : temperature (degrees Celsius)
        NO3 : nitrate concentration (umol/L)
        Si(OH)4 : silicon concentration (umol/L)
       
  Outputs: 
  
        printed data as follows : 
            1) Transferred Glucose Flux
            2) Cyanobacteria CO2 uptake
            3) Diatom CO2 uptake
            4) Diatom biomass generation
            5) Cyanobacteria biomass generation
            6) Carbon transfer value in table units (transferred C/basal cyano CO2)
            7) % of cyanobacterial metabolism supported by diatom
            8) % of basal cyanobacterial nitrogen fixation, will increase due to symbiosis

    This function will calculate the amount of nitrogen transferred by cyanobacteria and 
    carbon transferred by diatom given a certain temperature. This calculation is for
    comparison with Inomura et al.'s carbon cost analysis, to validate our treatment of
    symbiosis via FBA.

"""

def Nitrogen_CarbonTransfer(Inputs,T,NO3,Sili):

    import SymbiosisLoadModel, Extras
    
    [ML, DML, refmodel, IYSi, INSi, Inputs] = SymbiosisLoadModel.SymbiosisLoadModel(Inputs)

    # Determine henry's law coefficients for common atmospheric gases.
    Henry = Extras.henry(T)
    CO2sat = Henry[0]*(4.21e-4) # Partial pressure of CO2 in atm. 
    O2sat = Henry[1]*(.2105) # Partial pressure of O2 in atm. 
    N2sat = Henry[2]*(.781) # Partial pressure of N2 in atm. 

    
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

    # Temp-dependent Rubisco enzyme kinetics to set bounds of CO2 consumption for diatom.
    # Data can be found at:
    #   1) https://doi.org/10.1007/s11120-014-0067-8
    #   2) https://doi.org/10.1093/jxb/erw267
    ML.reactions.EX_co2_e.upper_bound = 1000
    ML.reactions.EX_co2_e.lower_bound = -Extras.RLPD(CO2sat*4,O2sat,T)
    ML.reactions.EX_hco3_e.lower_bound = 0

    # Michaelis-Menten saturation kinetics for nitrate uptake.
    # Data can be found at: https://doi.org/10.3354/meps09088
    ML.reactions.EX_no3_e.lower_bound = -((1.474*NO3)/(6.14+NO3))

    # Modified Michaelis-Menten saturation kinetics for silicic acid uptake, using a Hill coefficient
    # which represents cooperative allosteric interaction during silicic acid uptake. 
    # Data can be found at: https://doi.org/10.1104/pp.107.107094
    n = 1.9
    ML.reactions.EX_sio4h4_e.lower_bound = -((1.961*(Sili**n))/((8.1**n)+(Sili**n)))


    # Temp-dependent Rubisco enzyme kinetics to set bounds of CO2 consumption for cyanobacteria.
    # Data can be found at:
    #   1) https://doi.org/10.1007/s11120-014-0067-8
    #   2) https://doi.org/10.1093/jxb/erw267
    ML.reactions.get_by_id("EX_cpd00011[e]").upper_bound = 1000
    ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound = -Extras.RLPC(CO2sat*4,O2sat,T)
    ML.reactions.get_by_id("EX_cpd00011[e]").upper_bound = -Extras.RLPC(CO2sat*4,O2sat,T)*0.9999
    ML.reactions.get_by_id("EX_cpd00528[e]").lower_bound = -(0.9*N2sat)/(165+N2sat)
    ML.reactions.EX_nh4_e.lower_bound = 0
    ML.reactions.EX_nh4_e.upper_bound = 0

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
        ub=-10*refsolution.fluxes["EX_cpd00528[e]"])
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

    # Print all relevant values before ending function.
    print("Transferred Glucose: ",solution.fluxes["carbonT"])
    print("Cyanobacteria CO2 Intake: ",ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound)
    print("Diatom CO2 Intake: ",ML.reactions.EX_co2_e.lower_bound)
    print("Diatom Biomass Flux: ",solution.fluxes["DM_biomass_c"])
    print("Cyanobacteria Biomass Flux: ",solution.fluxes["biomass_eq_33047__vc"])
    print("Saved Carbon Transfer Value: ",(6*solution.fluxes["carbonT"]/(-ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound)))
    print("Cyanobacteria Metabolism %: ",100*((solution.fluxes["carbonT"]/((solution.fluxes["EX_cpd00011[e]"]/6)+solution.fluxes["carbonT"]))))
    print("Cyanobacteria Nitrogen %: ",100*((solution.fluxes["nh4T"]-refsolution.fluxes["EX_cpd00528[e]"])/(refsolution.fluxes["EX_cpd00528[e]"])))
    
    return