def MergeLoadModel(Inputs): # model: cyanobacterium Synechocystis sp. PCC 6803; model2: diatom Phaeodactylum tricornutum CCAP 1055/1
    import cobra
    import numpy as np
    from scipy.integrate import simpson
    from scipy.interpolate import RegularGridInterpolator
    import Light, Extras
    
    model = cobra.io.read_sbml_model('Models/Thaps_ML_model_n.xml')
    model2 = cobra.io.read_sbml_model('Models/iAnC892.xml')

    DML = cobra.io.read_sbml_model('Models/Thaps_ML_model_n.xml')
    refmodel = model2
    ML = model.merge(model2,objective='sum')
    ML.objective = {ML.reactions.get_by_id("DM_biomass_c"):1,ML.reactions.get_by_id("biomass_eq_33047__vc"):1}
    
    # Create reaction objects for each of the new reactions we are importing.
    # The first reaction is created implicitly below as an exchange flux for stored chryso.
    # The two here are a "transport" reaction, unstoring the chryso, and a conversion directly to glucose.
    nh4T = cobra.Reaction('nh4T')
    carbonT = cobra.Reaction('carbonT')
    G6PD_c = cobra.Reaction('G6PD_c')
    
    ML.add_reactions([nh4T,carbonT])
    ML.reactions.carbonT.build_reaction_from_string('glc__D_c --> cpd00027[vc]')
    ML.reactions.nh4T.build_reaction_from_string('cpd00013[hc] --> nh4_c')
    
    ML.add_reactions([G6PD_c])
    ML.reactions.G6PD_c.build_reaction_from_string('g6p_c --> glc__D_c')
    
    ML.reactions.nh4T.lower_bound = -1000
    ML.reactions.nh4T.upper_bound = 1000
    ML.reactions.carbonT.lower_bound = -1000
    ML.reactions.carbonT.upper_bound = 1000

    # Create reaction objects for each of the new reactions we are importing.
    # The first reaction is created implicitly below as an exchange flux for stored chryso.
    # The two here are a "transport" reaction, unstoring the chryso, and a conversion directly to glucose.
    gluconv = cobra.Reaction('gluconv')
    chrysoT = cobra.Reaction('chrysoT')
    GLCt1 = cobra.Reaction('GLCt1')
    
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

    ML.add_reactions([GLCt1,chrysoT,gluconv])
    ML.add_boundary(ML.metabolites.get_by_id("glc__D_e"), type="exchange")
    ML.add_boundary(model.metabolites.get_by_id("chryso_e"), type="exchange")
    ML.reactions.gluconv.build_reaction_from_string('chryso_c --> glc__D_c')
    ML.reactions.chrysoT.build_reaction_from_string('chryso_e --> chryso_c')
    ML.reactions.GLCt1.build_reaction_from_string('glc__D_e --> glc__D_c')
    
    DML.add_reactions([GLCt1,chrysoT,gluconv])
    DML.add_boundary(ML.metabolites.get_by_id("glc__D_e"), type="exchange")
    DML.add_boundary(model.metabolites.get_by_id("chryso_e"), type="exchange")
    DML.reactions.gluconv.build_reaction_from_string('chryso_c --> glc__D_c')
    DML.reactions.chrysoT.build_reaction_from_string('chryso_e --> chryso_c')
    DML.reactions.GLCt1.build_reaction_from_string('glc__D_e --> glc__D_c')
    
    # Remove author's biomass constraint forcing a growth rate at or above observed.
    # If this biomass constraint is used, all solutions at low light are infeasible.
    ML.reactions.biomass_eq_33047__hc.lower_bound = 0

    # Must run this cell or it will rely entirely on heterotrophic growth!!!
    # Defaults for co2 exchange are ub = 0, lb = 0.
    ML.reactions.EX_co2_e.upper_bound = 1000
    ML.reactions.EX_co2_e.lower_bound = -1000
    DML.reactions.EX_co2_e.upper_bound = 1000
    DML.reactions.EX_co2_e.lower_bound = -1000    

    # Must run this cell in order for photosynthesis to occur!
    # Defaults for light intensity are ub = 0, lb = 0. 
    # ML.reactions.EX_photon_e.lower_bound = -0.62
    # ML.reactions.EX_photon_e.upper_bound = 1000
    
    # Run this cell if you wish to perform autotrophic growth only.
    ML.reactions.EX_glc__D_e.lower_bound = 0
    ML.reactions.EX_glc__D_e.upper_bound = 0
    DML.reactions.EX_glc__D_e.lower_bound = 0
    DML.reactions.EX_glc__D_e.upper_bound = 0
    DML.reactions.EX_co2_e.upper_bound = 1000
    DML.reactions.EX_co2_e.lower_bound = -1000
    
    # Shut down or turn on the chrysolaminarin reaction.
    ML.reactions.EX_chryso_e.lower_bound = 0
    
    Iphoto = Inputs['Iphoto']
    
    # Weird photon input of this model.
    for k in Iphoto.keys():
        ML.reactions.get_by_id(k).upper_bound = 0
        ML.reactions.get_by_id(k).lower_bound = Iphoto[k]*-1.
        ML.reactions.get_by_id(k).upper_bound = Iphoto[k]*-0.9999
    
    ML.reactions.EX_no3_e.lower_bound = -0.001
    ML.reactions.EX_no2_e.lower_bound = 0
    ML.reactions.EX_nh4_e.lower_bound = 0
    ML.reactions.EX_pi_e.lower_bound = -1
    ML.reactions.EX_so4_e.lower_bound = -1
    ML.reactions.EX_sio4h4_e.lower_bound = -1
    ML.reactions.EX_cncbl3_e.lower_bound = -1
    ML.reactions.EX_co2_e.lower_bound = -10
    ML.reactions.EX_hco3_e.lower_bound = 0
    
    DML.reactions.EX_no3_e.lower_bound = -0.001
    DML.reactions.EX_no2_e.lower_bound = 0
    DML.reactions.EX_nh4_e.lower_bound = 0
    DML.reactions.EX_pi_e.lower_bound = -1
    DML.reactions.EX_so4_e.lower_bound = -1
    DML.reactions.EX_sio4h4_e.lower_bound = -1
    DML.reactions.EX_cncbl3_e.lower_bound = -1
    DML.reactions.EX_co2_e.lower_bound = -10
    DML.reactions.EX_hco3_e.lower_bound = 0    

    datamin_re = np.loadtxt("CarbonTransfer.csv",delimiter=',')
    Tn = np.linspace(-3,26,30)
    NO3 = np.linspace(0,.5,3)
    Sili = np.linspace(0.2,2,18)
    datamin = datamin_re.reshape(datamin_re.shape[0], datamin_re.shape[1] // len(Sili), len(Sili))
    Inputs['datamin'] = datamin
    IYSi = RegularGridInterpolator((Tn,NO3),datamin[:,:,-1])
    INSi = RegularGridInterpolator((Tn,NO3,Sili),datamin)
    
    return ML, DML, refmodel, IYSi, INSi, Inputs