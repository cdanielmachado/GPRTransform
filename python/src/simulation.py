from framed.analysis.simulation import pFBA, MOMA, lMOMA
from framed.analysis.variability import FVA
from framed.experimental.benchmark import evaluate, run_method

run_pFBA = lambda *args, **kwargs: run_method('pFBA', *args, **kwargs)
run_MOMA = lambda *args, **kwargs: run_method('MOMA', *args, **kwargs)
run_lMOMA = lambda *args, **kwargs: run_method('lMOMA', *args, **kwargs)


def gene_pFBA(model, objective=None, minimize=False, constraints=None, solver=None):
    reactions = [r_id for r_id in model.reactions if r_id.startswith('u_')]
    
    if constraints:
        constraints = model.convert_constraints(constraints)
    
    sol = pFBA(model, objective, minimize, constraints, reactions, solver)
    sol.values_converted = model.convert_fluxes(sol.values)
    return sol


def gene_MOMA(model, reference=None, constraints=None, solver=None):
    reactions = [r_id for r_id in model.reactions if r_id.startswith('u_')]
    
    if constraints:
        constraints = model.convert_constraints(constraints)
        
    if not reference:
        wt_sol = gene_pFBA(model)
        reference = wt_sol.values
        
    reference = {r_id: reference[r_id] for r_id in reactions}
    
    sol = MOMA(model, reference, constraints, solver)
    sol.values_converted = model.convert_fluxes(sol.values)
    return sol


def gene_lMOMA(model, reference=None, constraints=None, solver=None):
    reactions = [r_id for r_id in model.reactions if r_id.startswith('u_')]
    
    if constraints:
        constraints = model.convert_constraints(constraints)
    
    if not reference:
        wt_sol = gene_pFBA(model)
        reference = wt_sol.values
        
    reference = {r_id: reference[r_id] for r_id in reactions}
    
    sol = lMOMA(model, reference, constraints, solver)
    sol.values_converted = model.convert_fluxes(sol.values)
    return sol


def gene_FVA(model, obj_percentage=0, reactions=None, constraints=None):
    if constraints:
        constraints = model.convert_constraints(constraints)
        
    if not reactions:
        reactions = [r_id for r_id in model.reactions if r_id.startswith('u_')]

    result = FVA(model, obj_percentage, reactions, constraints)
    
    return result


def run_gpr_method(method, model, constraints, dataset, condition, **kwargs):

    genes = dataset.get_gene_deletions(condition)

    if not constraints:
        constraints = {}

    for gene in genes:
        if gene in model.genes:
            constraints['u_' + gene] = (0,0)

    if method in ['MOMA', 'lMOMA', 'ROOM']:
        reference = None
        ref_condition = kwargs.get('reference_condition', None)
        if ref_condition:
            run_pFBA = lambda *args, **kwargs: run_gpr_method('pFBA', *args, **kwargs)
            _, _, sol = evaluate(run_pFBA, model, dataset, ref_condition, **kwargs)
            reference = sol.values
        else:
            print 'Warning: Must specify reference condition for MOMA/lMOMA/ROOM simulations.'

    if method == 'pFBA':
        sol = gene_pFBA(model, constraints=constraints)
    elif method == 'MOMA':
        sol = gene_MOMA(model, constraints=constraints, reference=reference)
    elif method == 'lMOMA':
        sol = gene_lMOMA(model, constraints=constraints, reference=reference)
    else:
        sol = None

    return sol, model.convert_fluxes(sol.values)


run_g_pFBA = lambda *args, **kwargs: run_gpr_method('pFBA', *args, **kwargs)
run_g_MOMA = lambda *args, **kwargs: run_gpr_method('MOMA', *args, **kwargs)
run_g_lMOMA = lambda *args, **kwargs: run_gpr_method('lMOMA', *args, **kwargs)


