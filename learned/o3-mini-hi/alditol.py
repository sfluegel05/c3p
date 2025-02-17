"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: alditol
Alditols are acyclic polyols having the general formula HOCH2[CH(OH)]nCH2OH 
(formally derivable from an aldose by reduction of the carbonyl group).
This script looks for a contiguous carbon chain where the end carbons
appear as primary alcohols (CH2OH) and the interior carbons as secondary alcohols (CH(OH)).
"""
from rdkit import Chem
from collections import deque

def get_connected_components(graph):
    """
    Given a graph as a dictionary {node: [neighbors]}, returns a list of sets,
    each being a connected component.
    """
    seen = set()
    components = []
    for node in graph:
        if node not in seen:
            comp = set()
            queue = deque([node])
            seen.add(node)
            while queue:
                current = queue.popleft()
                comp.add(current)
                for nbr in graph[current]:
                    if nbr not in seen:
                        seen.add(nbr)
                        queue.append(nbr)
            components.append(comp)
    return components

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is defined as an acyclic polyol with the general pattern
    HOCH2[CH(OH)]nCH2OH (reduction product of an aldose).
    
    This routine first adds explicit hydrogens so that -OH groups are detectable.
    It then finds carbons that bear at least one hydroxyl group. From these,
    it looks for a connected subgraph (chain) that is linear: two terminal carbons
    (primary alcohols, expected to have two bound hydrogens) and interior carbons
    (secondary alcohols, expected to have one hydrogen).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as alditol, False otherwise.
        str : Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to count them in alcohol groups
    mol = Chem.AddHs(mol)

    candidate_atoms = []  # indices of carbon atoms that bear at least one -OH group
    # For each carbon atom, check if it has a neighbor oxygen that is bound to at least one hydrogen.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # only consider carbons
        oh_found = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # check if oxygen has at least one hydrogen neighbor.
                h_count = sum(1 for nn in nbr.GetNeighbors() if nn.GetAtomicNum() == 1)
                if h_count >= 1:
                    oh_found = True
                    break
        if oh_found:
            candidate_atoms.append(atom.GetIdx())
    
    if not candidate_atoms:
        return False, "No carbon atoms with alcohol (-OH) groups detected"

    # Build a graph among candidate carbons: nodes are candidate atom indices, edges exist if a bond exists between them.
    candidate_set = set(candidate_atoms)
    graph = {idx: [] for idx in candidate_set}
    for idx in candidate_set:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in candidate_set:
                graph[idx].append(nbr_idx)
    
    # Look for connected components in this candidate graph.
    components = get_connected_components(graph)
    # We search for a component that is at least 2 carbons long (to have terminal groups)
    for comp in components:
        if len(comp) < 2:
            continue
        # In a linear (chain) structure, exactly two atoms are terminal (i.e. degree 1)
        degs = [len(graph[node]) for node in comp]
        terminal_count = sum(1 for d in degs if d == 1)
        # The ideal chain: two terminal (degree-1) nodes and the remainder are degree 2.
        if terminal_count != 2:
            continue
        # Now check that the terminal carbons appear as primary alcohols:
        # Here primary means that the carbon likely is CH2–OH (we expect two hydrogens on it).
        terminal_nodes = [node for node in comp if len(graph[node]) == 1]
        terminals_ok = True
        for node in terminal_nodes:
            atom = mol.GetAtomWithIdx(node)
            h_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
            # Expect at least 2 hydrogens for a CH2 group.
            if h_count < 2:
                terminals_ok = False
                break
        if not terminals_ok:
            continue
        # Check interior nodes (if any) are secondary alcohols: expect exactly one hydrogen attached.
        interior_nodes = [node for node in comp if len(graph[node]) == 2]
        interiors_ok = True
        for node in interior_nodes:
            atom = mol.GetAtomWithIdx(node)
            h_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
            # Often a secondary carbon (CH) has one hydrogen.
            if h_count < 1:
                interiors_ok = False
                break
        if not interiors_ok:
            continue

        # Optional: one may also require that the chain is acyclic.
        # Even if the whole molecule may contain rings elsewhere, here our candidate chain should be an open chain.
        # We also check that none of the candidate atoms is in a ring.
        in_ring = any(mol.GetAtomWithIdx(node).IsInRing() for node in comp)
        if in_ring:
            continue

        reason = f"Found alditol chain with {len(comp)} contiguous carbon atoms carrying -OH groups (terminal primary and interior secondary alcohols)."
        return True, reason

    return False, "No linear alditol chain pattern detected"

# Example usage (uncomment to run tests):
# test_smiles = "OC[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)CO"  # L-iditol example
# result, msg = is_alditol(test_smiles)
# print(result, msg)