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
with open('zinc.pickle', 'rb') as z:
    zinc = load(z)


def b64(m):
    return encodebytes(m.depict().encode()).decode().replace('\n', '')


def visualization(G, target):
    template = 'data:image/svg+xml;base64,'

    synky = list(G._succ[target])
    if synky:
        G.remove_nodes_from(synky)

    info('go to visualization')
    data = {}
    _all = {target: str(target)}
    nodes = [{'id': str(target), 'image': template + b64(target), 'shape': 'image', 'size': 80,
              "shapeProperties": {"useBorderWithImage": True}, "color": {"border": "red", "background": "white"}}]
    with db_session:
        for node in G.nodes():
            if node != target:
                if isinstance(node, int):
                    description = G.node[node]['data']
                    _all.update({node: description})
                    item = {'id': description, 'label': 'reaction',
                            'color': {'border': 'black',
                                      'background': 'white'},
                            'shape': 'box',
                            'widthConstraint': {'maximum': 300}, 'font': {'size': 20}}
                else:
                    _all.update({node: str(node)})
                    color_border = "#339FFF"
                    if 'zinc' in G.nodes[node]:
                        color_border = "green"
                    item = {'id': str(node), 'image': template + b64(node), 'shape': 'image', 'size': 70,
                            "shapeProperties": {"useBorderWithImage": True},
                            "color": {"border": color_border, "background": "white"}}
                nodes.append(item)
        data['nodes'] = nodes
        edges = []
        for source, target in G.edges():
            item = {'width': 2, 'from': _all[source], 'to': _all[target]}
            if isinstance(source, int):
                item['color'] = {'color': 'red'}
            else:
                item['color'] = {'color': 'green'}
            edges.append(item)
        data['edges'] = edges
        info('Finish Him!')
    return dumps(data)


def similarity(target):
    molecules = list()
    similar_molecules = Molecule.find_similar(target, page=1, pagesize=100)
    for pair in similar_molecules:
        molecule, tan = pair
        if MoleculeReaction.exists(molecule=molecule, is_product=True) and len(molecules) < 10:
            molecules.append((molecule, round(tan, 2), 0))

    if not molecules:
        info('Molecule-target and similar molecules not found as product')
        raise Exception

    return molecules


def substructure(target):
    molecules = list()
    rctns = None
    n = 0
    _target = Molecule.find_structure(target)
    if _target:
        if MoleculeReaction.exists(molecule=_target, is_product=True):
            reactions = [x.id for x in _target.reactions_entities(pagesize=100, product=True)]
            shuffle(reactions)
            molecules.append((_target, 1, reactions[:10]))

    substructures = Molecule._structure_query(target, 'substructure')
    for pair in substructures:
        if n and rctns and len(molecules) < 10 and m:
            shuffle(rctns)
            molecules.append((m, round(tan, 2), rctns[:10]))
        m = None
        n += 1
        sub, tan = pair
        if sub != _target:
            for mapping in target.get_substructure_mapping(sub.structure, limit=0):
                if MoleculeReaction.exists(molecule=sub, is_product=True):
                    reactions = sub.reactions_entities(pagesize=100, product=True)
                    rctns = []
                    for r in reactions:
                        for mppng in select(x for x in MoleculeReaction.entity if x.molecule == sub and x.reaction == r):
                            mppng = mppng.mapping
                            if not set(r.cgr.center_atoms).isdisjoint(mppng.get(x, x) for x in mapping.values()):
                                m = sub
                                rctns.append(r.id)

    if molecules:
        return molecules
    else:
        info('Substructures not found as product')
        raise Exception


def update_graph(target, seen, stack, available, g, n=0, added_reactions=None):
    if added_reactions is None:
        added_reactions = set()

    # create graph
    while stack:
        mt, st = stack.pop(0)
        st -= 1
        reactions = mt.reactions_entities(pagesize=100, product=True)
        shuffle(reactions)
        for r in reactions[:10]:
            if r.id in added_reactions:
                continue
            data = list(r.metadata)[0].data
            reagents = ''
            if data['reagents']:
                reagents = f"Reagents: {(data['reagents'])}. "
            g.add_node(n, data=f"{data['source_id']}; {data['text']}".replace('+\n', '') +
                               f" {reagents}Registration year: {data['year']}")
            added_reactions.add(r.id)
            for m in r.molecules:
                obj_mol = m.molecule
                structure = obj_mol.structure
                atoms = [atom[-1].element for atom in structure.atoms()]
                if not m.is_product:
                    g.add_edge(structure, n)
                    if st and obj_mol not in seen:
                        seen.add(obj_mol)
                        if available:
                            if mt == target and bytes(structure) not in zinc \
                                    and not (atoms.count('C') <= 2 or len(atoms) <= 6):
                                stack.append((obj_mol, st))
                            elif bytes(mt.structure) not in zinc and bytes(structure) not in zinc \
                                    and not (atoms.count('C') <= 2 or len(atoms) <= 6):
                                stack.append((obj_mol, st))
                        else:
                            stack.append((obj_mol, st))
                else:
                    g.add_edge(n, structure)
                if bytes(structure) in zinc or atoms.count('C') <= 2 or len(atoms) <= 6:
                    g.nodes[structure]['zinc'] = 1
            n += 1


def paths_of_synthesis_for_target_molecule(target, available, stages):
    g = DiGraph()
    update_graph(target, {target}, [(target, stages)], available, g)
    if not len(g):

        raise Exception

    return visualization(g, target.structure)


def paths_of_synthesis_for_substructure_search(target, available, reactions, stages):
    seen = set()
    n = 0
    g = DiGraph()
    stack = []
    added_reactions = set()
    for id in reactions:
        reaction = Reaction[id]
        if id not in added_reactions:
            data = list(reaction.metadata)[0].data
            reagents = ''
            if data['reagents']:
                reagents = f"Reagents: {(data['reagents'])}; "
            g.add_node(n, data=f"{data['source_id']}; {data['text']}".replace('+\n', '') +
                               f" {reagents}registration year: {data['year']}")
            added_reactions.add(id)
        else:
            continue
        for m in reaction.molecules:
            obj_mol = m.molecule
            structure = obj_mol.structure
            atoms = [atom[-1].element for atom in structure.atoms()]
            if not m.is_product:
                if stages and obj_mol not in seen:
                    seen.add(obj_mol)
                    if available:
                        if bytes(structure) not in zinc and not (atoms.count('C') <= 2 or len(atoms) <= 6):
                            stack.append((obj_mol, stages))
                    else:
                        stack.append((obj_mol, stages))
                g.add_edge(structure, n)
            else:
                g.add_edge(n, structure)
            if bytes(structure) in zinc or atoms.count('C') <= 2 or len(atoms) <= 6:
                g.nodes[structure]['zinc'] = 1
        n += 1

    update_graph(target, seen, stack, available, g, n, added_reactions)
    target = target.structure
    return visualization(g, target)
