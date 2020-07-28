from base64 import encodebytes
from CGRdb import Molecule, Reaction
from CGRdb.database import MoleculeReaction
from json import dumps
from logging import info, basicConfig, INFO
from networkx import DiGraph
from pickle import load
from pony.orm import db_session, select
from random import shuffle

basicConfig(level=INFO)
with open('zinc_reaxys.pickle', 'rb') as z:
    zinc = load(z)


def b64(m):
    return encodebytes(m.depict().encode()).decode().replace('\n', '')


def visualization(g, target):
    template = 'data:image/svg+xml;base64,'
    info('go to visualization')
    data = {}
    _all = {target: str(target)}
    nodes = [{'id': str(target), 'image': template + b64(target), 'shape': 'image', 'size': 80,
              "shapeProperties": {"useBorderWithImage": True}, "color": {"border": "red", "background": "white"}}]
    for node in g.nodes():
        if node != target:
            if isinstance(node, int):
                description = g.nodes[node]['data']
                _all.update({node: description})
                item = {'id': description, 'label': 'reaction',
                        'color': {'border': 'black',
                                  'background': 'white'},
                        'shape': 'box',
                        'widthConstraint': {'maximum': 300}, 'font': {'size': 20}}
            else:
                _all.update({node: str(node)})
                color_border = "#339FFF"
                if 'zinc' in g.nodes[node]:
                    color_border = "green"
                item = {'id': str(node), 'image': template + b64(node), 'shape': 'image', 'size': 70,
                        "shapeProperties": {"useBorderWithImage": True},
                        "color": {"border": color_border, "background": "white"}}
            nodes.append(item)
    data['nodes'] = nodes
    edges = []
    for source, target in g.edges():
        item = {'width': 2, 'from': _all[source], 'to': _all[target]}
        if isinstance(source, int):
            item['color'] = {'color': 'red'}
        else:
            item['color'] = {'color': 'green'}
        edges.append(item)
    data['edges'] = edges
    info('Finish Him!')
    return dumps(data)


def similarity(target, cutoff=10):
    molecules = list()
    similar_molecules = Molecule.find_similar(target, page=1, pagesize=100)
    for pair in similar_molecules:
        molecule, tan = pair
        if MoleculeReaction.exists(molecule=molecule, is_product=True) and len(molecules) < cutoff:
            molecules.append((molecule, round(tan, 2), 0))

    if not molecules:
        raise Exception('Molecule-target and similar molecules not found as product')

    return molecules


def substructure(target, cutoff=10):
    molecules = list()

    _target = Molecule.find_structure(target)
    if _target and MoleculeReaction.exists(molecule=_target, is_product=True):
        reactions = [x.id for x in _target.reactions_entities(pagesize=100, product=True)]
        shuffle(reactions)
        molecules.append((_target, 1, reactions[:cutoff]))

    n = 0
    rctns = None
    substructures = Molecule.find_substructures(target)
    if substructures:
        for sub, tan in zip(substructures.molecules(), substructures.tanimotos()):
            if n and rctns and len(molecules) < cutoff and mol:
                shuffle(rctns)
                molecules.append((mol, round(tan, 2), rctns[:cutoff]))
            mol = None
            n += 1
            if sub != _target:
                for mapping in target.get_mcs_mapping(sub.structure, limit=0):
                    if MoleculeReaction.exists(molecule=sub, is_product=True):
                        reactions = sub.reactions_entities(pagesize=100, product=True)
                        rctns = []
                        for r in reactions:
                            for mppng in select(x for x in MoleculeReaction.entity if x.molecule == sub and
                                                                                      x.reaction == r):
                                mppng = mppng.mapping
                                if not set(r.cgr.center_atoms).isdisjoint(mppng.get(x, x) for x in mapping.values()):
                                    mol = sub
                                    rctns.append(r.id)

    if molecules:
        return molecules
    else:
        raise Exception('Substructures not found as product')


def update_graph(target, seen, stack, available, g, cutoff, n=0, added_reactions=None):
    if added_reactions is None:
        added_reactions = set()

    # create graph
    while stack:
        mt, st = stack.pop(0)
        st -= 1
        reactions = mt.reactions_entities(pagesize=100, product=True)
        shuffle(reactions)
        for r in reactions[:cutoff]:
            if r.id in added_reactions:
                continue
            data = list(r.metadata)[0].data
            added_reactions.add(r.id)
            products = set()
            reactants = set()
            for m in r._molecules:
                if not m.is_product:
                    reactants.add(m.molecule)
                else:
                    products.add(m.molecule)

            if target in reactants:
                continue
            else:
                rem = reactants.intersection(products)
                for obj_mol in reactants:
                    if obj_mol in rem:
                        continue
                    structure = obj_mol.structure
                    atoms = [atom[-1].atomic_symbol for atom in structure.atoms()]
                    g.add_edge(structure, n)
                    if st and obj_mol not in seen:
                        seen.add(obj_mol)
                        if available:
                            if bytes(structure) not in zinc and not (atoms.count('C') <= 2 or len(atoms) <= 6):
                                if mt == target:
                                    stack.append((obj_mol, st))
                                elif bytes(mt.structure) not in zinc:
                                    stack.append((obj_mol, st))
                            else:
                                g.nodes[structure]['zinc'] = 1
                        else:
                            stack.append((obj_mol, st))
                            if bytes(structure) in zinc and (atoms.count('C') <= 2 or len(atoms) <= 6):
                                g.nodes[structure]['zinc'] = 1
                for obj_mol in products:
                    if obj_mol in rem:
                        continue
                    g.add_edge(n, obj_mol.structure)
                g.nodes[n]['data'] = f"{data}"
                n += 1


def paths_of_synthesis_for_target_molecule(target, available, stages, cutoff=10):
    g = DiGraph()
    update_graph(target, {target}, [(target, stages)], available, g, cutoff)
    if not len(g):

        raise Exception

    return visualization(g, target.structure)


def paths_of_synthesis_for_substructure_search(target, available, reactions, stages, cutoff=10):
    seen = set()
    n = 0
    g = DiGraph()
    stack = []
    added_reactions = set()
    for idd in reactions:
        reaction = Reaction[idd]
        if idd in added_reactions:
            continue
        data = list(reaction.metadata)[0].data
        added_reactions.add(idd)

        reactants = set()
        products = set()
        for m in reaction._molecules:
            if not m.is_product:
                reactants.add(m.molecule)
            else:
                products.add(m.molecule)

        if target in reactants:
            continue
        else:
            rem = reactants.intersection(products)
            for obj_mol in reactants:
                if obj_mol in rem:
                    continue
                structure = obj_mol.structure
                atoms = [atom[-1].atomic_symbol for atom in structure.atoms()]
                g.add_edge(structure, n)
                if stages and obj_mol not in seen:
                    seen.add(obj_mol)
                    if available:
                        if bytes(structure) not in zinc and not (atoms.count('C') <= 2 or len(atoms) <= 6):
                            stack.append((obj_mol, stages))
                        else:
                            g.nodes[structure]['zinc'] = 1
                    else:
                        stack.append((obj_mol, stages))
                        if bytes(structure) in zinc and (atoms.count('C') <= 2 or len(atoms) <= 6):
                            g.nodes[structure]['zinc'] = 1
            for obj_mol in products:
                if obj_mol in rem:
                    continue
                g.add_edge(n, obj_mol.structure)
            g.nodes[n]['data'] = f"{data}"

    update_graph(target, seen, stack, available, g, cutoff, n=n, added_reactions=added_reactions)
    target = target.structure
    return visualization(g, target)
