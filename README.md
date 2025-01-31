# GCB_Diatom
Companion code for submitted GCB article (DOI to come).

This code represents a novel application of dynamic flux balance analysis (dFBA) for metabolic modeling of marine diatoms and cyanobacteria. 
The biochemical kinetics of these organisms are very important to the global carbon cycle, as diatoms are responsible for up to a quarter of photosynthesis (thus, carbon sequestration) and cyanobacteria can help to relieve diatom's nitrogen stress.

We apply a hierarchical optimization approach, originally proposed in Hoffner, 2013. This approach guarantees that any fluxes calculated within FBA are unique, and will not lead to convergence/stability issues for the dynamic differential equation solver.

In addition, we have ported the Simple Model of Atmospheric Radiation and Transfer of Sunshine (SMARTS2) from Gueymard, 1995 from FORTRAN to Python for our use.

Our model incorporates temperature-dependent RuBisCo-limited kinetics, as well as correlations for all relevant limiting nutrients (nitrogen, silicic acid, iron). 
Our peak biomass concentrations agree with known biomass concentrations during diatom blooms. 
Our results suggest that climate change may lead to a positive feedback loop, in which marine phytoplankon consume less carbon with increasing sea temperatures.
