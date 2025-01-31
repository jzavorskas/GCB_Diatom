def MergeTransfer(Inputs,y,T,ML,IYSi,INSi):
    
    import numpy as np
    import cobra
    
    Tn = np.linspace(-3,26,30)
    
    if y[6] > 0.75:
        CT = np.interp(T,Tn,Inputs['datamin'][:,0,-1])
        CarbonConstraint = ML.problem.Constraint(
            ML.reactions.carbonT.flux_expression,
            lb=-(CT/6)*ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound,
            ub=-(CT/6)*ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound)
    elif y[6] > 0.32:
        CT = INSi([T,0,y[5]])
        CarbonConstraint = ML.problem.Constraint(
            ML.reactions.carbonT.flux_expression,
            lb=-(CT/6)*ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound,
            ub=-(CT/6)*ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound)
    else:
        CarbonConstraint = ML.problem.Constraint(
            ML.reactions.carbonT.flux_expression,
            lb=0,
            ub=0)
        CT = 0
        
    return CarbonConstraint, -(CT/6)*ML.reactions.get_by_id("EX_cpd00011[e]").lower_bound