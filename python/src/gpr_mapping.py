from framed.core.model import Metabolite, Reaction, Compartment
from framed.core.cbmodel import GPRAssociation
from framed.core.transformation import make_irreversible


def split_isozymes(model):

    mapping = dict()

    for r_id, reaction in model.reactions.items():
                
        gpr = model.gpr_associations[r_id]
        
        if gpr is not None and len(gpr.proteins) > 1:
            mapping[r_id] = []
            for i, protein in enumerate(gpr.proteins): 
                r_id_new = '{}_pc{}'.format(reaction.id, i+1)
                mapping[r_id].append(r_id_new)
                reaction_new = Reaction(r_id_new, reaction.name, reaction.reversible, reaction.stoichiometry, reaction.regulators)
                model.add_reaction(reaction_new)
                model.set_flux_bounds(r_id_new, *model.bounds[r_id])
                model.set_reaction_objective(r_id_new, model.objective[r_id])
                gpr_new = GPRAssociation()
                gpr_new.proteins.append(protein)
                model.set_gpr_association(r_id_new, gpr_new)
            model.remove_reaction(r_id)
    return mapping


def genes_to_species(model):
    spontaneous = 's0001'
    new_reactions = []
    compartment = Compartment('genes', 'gene pool')
    model.add_compartment(compartment)
    
    for gene in model.genes.values():
        if gene.id == spontaneous:
            continue
        model.add_metabolite(Metabolite(gene.id, gene.id, 'genes'))
        r_id = 'u_{}'.format(gene.id)
        reaction = Reaction(r_id, r_id, False, {gene.id: 1})
        model.add_reaction(reaction)
        new_reactions.append(r_id)
    
    for r_id, reaction in model.reactions.items():
        gpr = model.gpr_associations[r_id]

        if gpr is not None:
            if len(gpr.proteins) > 1:
                print 'error: isozymes not split:', r_id
                return
            elif len(gpr.proteins) == 1:
                for g_id in gpr.proteins[0].genes:
                    if g_id != spontaneous:
                        reaction.stoichiometry[g_id] = -1
    
    return new_reactions     


def merge_fluxes(fluxes, mapping_rev, mapping_iso, net=True):
    fluxes = fluxes.copy()
    
    for r_id, r_ids in mapping_iso.items():
        fluxes[r_id] = sum([fluxes[r_id2] for r_id2 in r_ids])
        for r_id2 in r_ids:
            del fluxes[r_id2]
    
    for r_id, (fwd_id, bwd_id) in mapping_rev.items():
        if net:
            fluxes[r_id] = fluxes[fwd_id] - fluxes[bwd_id]
        else:
            fluxes[r_id] = fluxes[fwd_id] + fluxes[bwd_id]
        del fluxes[fwd_id]
        del fluxes[bwd_id]
        
    return fluxes


def convert_constraints(constraints, mapping_rev, mapping_iso):
    constraints = constraints.copy()
    
    for r_id, (fwd_id, bwd_id) in mapping_rev.items():
        if r_id in constraints:
            x = constraints[r_id]
            lb, ub = x if isinstance(x, tuple) else (x, x)            
            lb_fwd = max(0, lb) if lb is not None else 0
            ub_fwd = max(0, ub) if ub is not None else None
            lb_bwd = max(-ub, 0) if ub is not None else 0
            ub_bwd = max(-lb, 0) if lb is not None else None
            constraints[fwd_id] = (lb_fwd, ub_fwd)
            constraints[bwd_id] = (lb_bwd, ub_bwd)
            del constraints[r_id]

    for r_id, r_ids in mapping_iso.items():
        if r_id in constraints:
            x = constraints[r_id]
            ub = x[1] if isinstance(x, tuple) else x  
            for r_id2 in r_ids:
                constraints[r_id2] = (0, ub)
            del constraints[r_id]
    
    return constraints
            

def transform(model, inplace=True):
    if not inplace:
        model = model.copy()

    mapping_rev = make_irreversible(model)
    mapping_iso = split_isozymes(model)
    u_reactions = genes_to_species(model)
    model.convert_fluxes = lambda x, net=True: merge_fluxes(x, mapping_rev, mapping_iso, net)
    model.convert_constraints = lambda x: convert_constraints(x, mapping_rev, mapping_iso)
    return model


