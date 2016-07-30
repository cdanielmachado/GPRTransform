
from sympy.parsing.sympy_parser import parse_expr
from sympy import to_dnf, Or, And, Not
from framed.analysis.deletion import reaction_deletion
from framed.solvers.solver import Status
from framed.solvers import solver_instance
from multiprocessing import Pool, cpu_count
from pandas import DataFrame


class Callable:
    
    def __init__(self, func, state):
        self.func = func
        self.state = state
        
    def __call__(self, x):
        return self.func(x, self.state)


def parse_model_rules(model):
    parsed_rules = {}
    for r_id, rule in model.rules.items():
        if rule:
            rule = '( ' + rule.replace(' and ', ' & ').replace(' or ', ' | ').strip() + ' )'
            parsed_rules[r_id] = parse_expr(rule)
        else:
            parsed_rules[r_id] = None
    model.parsed_rules = parsed_rules



def reaction_to_gene_deletions(reaction_set, rules):    
    rule_set = [rules[r_id] for r_id in reaction_set]
    merged_rule = to_dnf(Not(Or(*rule_set)))
    
    gene_sets = []
    if type(merged_rule) is Or: 
        for sub_expr in merged_rule.args:
            gene_set = []
            if type(sub_expr) is And:
                gene_set = [str(not_gene.args[0]) if type(not_gene) is Not else 'error' 
                            for not_gene in sub_expr.args]
            else:
                gene_set = [str(sub_expr.args[0]) if type(sub_expr) is Not else 'error']
            gene_set = tuple(sorted(set(gene_set)))
            gene_sets.append(gene_set)
    elif type(merged_rule) is And:
        gene_set = [str(not_gene.args[0]) if type(not_gene) is Not else 'error' 
                    for not_gene in merged_rule.args]
        gene_set = tuple(sorted(set(gene_set)))
        gene_sets.append(gene_set)
    else:
        gene_set = [str(merged_rule.args[0]) if type(merged_rule) is Not else 'error']
        gene_set = tuple(sorted(set(gene_set)))
        gene_sets.append(gene_set)            
    return gene_sets


def expand_gene_sets(model, reaction_sets):
    
    rxns2genes = {}
    all_rxns = set()
    
    not_knockable = set([r_id for r_id, rule in model.rules.items() 
                         if rule == '' or 's0001' in rule])
    
    knockable_sets = []
    for reaction_set in reaction_sets:
        reaction_set = set(reaction_set)
        if len(reaction_set & not_knockable) == 0:
            all_rxns.update(reaction_set)
            knockable_sets.append(tuple(sorted(reaction_set)))
            
    my_func = Callable(reaction_to_gene_deletions, model.parsed_rules)
    p = Pool(cpu_count())
    gene_sets = p.map(my_func, knockable_sets)
    
    rxns2genes = dict(zip(knockable_sets, gene_sets))
    return rxns2genes


def genes_to_reactions(gene_set, model):
    
    genes_state = {gene: gene not in gene_set for gene in model.genes}
    del_rxns = [r_id for r_id, f in model.rule_functions.items() if not f(genes_state)]
                
    return tuple(sorted(del_rxns))


def gene_to_reaction_sets(model, rxns2genes): 
    rxns2rxns = {}
    for reaction_set, gene_sets in rxns2genes.items():
        my_func = Callable(genes_to_reactions, model)
        #p = Pool(cpu_count())
        rxn_sets = map(my_func, gene_sets)
        rxns2rxns[reaction_set] = set(rxn_sets)    
    
    return rxns2rxns


def gene_to_reaction_sets2(model, gene_sets): 
    gene2rxns = {}
    for gene_set in gene_sets:
        gene_set = tuple(sorted(gene_set))
        gene2rxns[gene_set] = genes_to_reactions(gene_set, model)
        
    return gene2rxns


def validate_reaction_deletions(model, del_rxns, constraints, solver=None):
    
    # if pFBA is too slow, replace it with FBA w/objective tilting !
    sol = reaction_deletion(model, del_rxns, method='pFBA', solver=solver)
    
    if sol.status == Status.OPTIMAL:
        sol.valid = all([lb < sol.values[r_id] < ub for r_id, (lb, ub) in constraints.items()])
    else:
        sol.valid = False
    return sol


def build_reaction_solution_pool(model, rxns2rxns, constraints=None):
    
    reaction_sets = reduce(set.__or__, rxns2rxns.values())
    
    solution_pool = {}
    solver = solver_instance()
    solver.build_problem(model)
    
    for reaction_set in reaction_sets:
        solution_pool[reaction_set] = validate_reaction_deletions(model, reaction_set, constraints, solver)
    
    return solution_pool


def build_reaction_solution_pool2(model, gene2rxns, constraints=None):
        
    solution_pool = {}
    solver = solver_instance()
    solver.build_problem(model)
    
    reaction_sets = set(gene2rxns.values())
    
    for reaction_set in reaction_sets:
        solution_pool[reaction_set] = validate_reaction_deletions(model, reaction_set, constraints, solver)
    
    return solution_pool
