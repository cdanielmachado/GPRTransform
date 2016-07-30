from __future__ import division
from framed.analysis.simulation import FBA, pFBA
from framed.experimental.benchmark import run_method
from framed.omics.simulation import GIMME, eflux
from numpy import percentile


def gpr_GIMME(model, gene_exp, cutoff=25, growth_frac=0.9, constraints=None, parsimonious=False):

    threshold = percentile(gene_exp.values(), cutoff)
    coeffs = {'u_' + gene: threshold-val for gene, val in gene_exp.items() if val < threshold}

    if constraints:
        constraints = model.convert_constraints(constraints)
    else:
        constraints = {}

    pre_solution = FBA(model, constraints=constraints)
    biomass = model.detect_biomass_reaction()
    constraints[biomass] = (growth_frac * pre_solution.values[biomass], None)

    if parsimonious:
        solution = pFBA(model, objective=coeffs, minimize=True, constraints=constraints)
    else:
        solution = FBA(model, objective=coeffs, minimize=True, constraints=constraints)

    solution.values_converted = model.convert_fluxes(solution.values)
    return solution


def gpr_eflux(model, gene_exp, scale_rxn, scale_value, constraints=None, parsimonious=False):

    max_exp = max(gene_exp.values())
    bounds = {}

    for r_id, (lb, ub) in model.bounds.items():
        if r_id.startswith('u_'):
            gene = r_id[len('u_'):]
            val = gene_exp[gene] / max_exp if gene in gene_exp else 1
            bounds[r_id] = (0, val)
        else:
            ub2 = 1 if ub is None or ub > 0 else 0
            bounds[r_id] = (0, ub2)

    if constraints:
        constraints = model.convert_constraints(constraints)
        for r_id, x in constraints.items():
            ub = x[1] if isinstance(x, tuple) else x 
            ub2 = 1 if ub is None or ub > 0 else 0
            bounds[r_id] = (0, ub2)

    if parsimonious:
        sol = pFBA(model, constraints=bounds)
    else:
        sol = FBA(model, constraints=bounds)

    sol.values_converted = model.convert_fluxes(sol.values)

    k = abs(scale_value / sol.values_converted[scale_rxn])

    for r_id, val in sol.values_converted.items():
        sol.values_converted[r_id] = val * k

    return sol



def run_GIMME(model, constraints, dataset, condition, **kwargs):
    growth_frac = kwargs.get('growth_frac', 0.9)
    cutoff = kwargs.get('cutoff', 25)
    parsimonious = kwargs.get('parsimonious', False)

    gene_exp = dataset.get_transcriptomics(conditions=condition)
    gene_exp = {gene: val for gene, val in gene_exp.iteritems() if gene in model.genes}
    sol = GIMME(model, gene_exp, growth_frac=growth_frac, cutoff=cutoff,
                constraints=constraints, parsimonious=parsimonious)
    return sol, sol.values


def run_gpr_GIMME(model, constraints, dataset, condition, **kwargs):
    growth_frac = kwargs.get('growth_frac', 0.9)
    cutoff = kwargs.get('cutoff', 25)
    parsimonious = kwargs.get('parsimonious', False)

    gene_exp = dataset.get_transcriptomics(conditions=condition)
    gene_exp = {gene: val for gene, val in gene_exp.iteritems() if gene in model.genes}

    sol = gpr_GIMME(model, gene_exp, growth_frac=growth_frac, cutoff=cutoff,
                    constraints=constraints, parsimonious=parsimonious)
    return sol, sol.values_converted


def run_eflux(model, constraints, dataset, condition, **kwargs):

    parsimonious = kwargs.get('parsimonious', False)

    gene_exp = dataset.get_transcriptomics(conditions=condition)
    gene_exp = {gene: val for gene, val in gene_exp.iteritems() if gene in model.genes}
    
    scale_rxn, scale_value = constraints.items()[0]

    if len(constraints) > 1:
        print 'warning: ambiguous scale reaction, using:', scale_rxn
    
    sol = eflux(model, gene_exp, scale_rxn, scale_value, constraints, parsimonious)
    return sol, sol.values


def run_gpr_eflux(model, constraints, dataset, condition, **kwargs):

    parsimonious = kwargs.get('parsimonious', False)

    gene_exp = dataset.get_transcriptomics(conditions=condition)
    gene_exp = {gene: val for gene, val in gene_exp.iteritems() if gene in model.genes}

    if len(constraints) > 1:
        print 'warning: ambiguous scale reaction'
    
    scale_rxn, scale_value = constraints.items()[0]
    
    sol = gpr_eflux(model, gene_exp, scale_rxn, scale_value, constraints, parsimonious)
    return sol, sol.values_converted


