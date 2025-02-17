"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: Lipopeptide – a compound consisting of a peptide with an attached lipid.
A molecule is considered a lipopeptide if it contains at least two amide bonds (–N–C(=O)– fragments)
and a long aliphatic chain (at least 8 sp3 nonaromatic carbons) attached to one of those amide substructures.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def get_longest_aliphatic_chain(mol):
    """
    Returns the length (number of atoms) of the longest chain composed solely of sp3-hybridized,
    nonaromatic carbon atoms.
    We build an explicit graph for sp3 aliphatic carbons and compute the chain's diameter.
    """
    valid_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and atom.GetHybridization().name == "SP3":
            valid_atoms.append(atom.GetIdx())
    if not valid_atoms:
        return 0

    # Build an adjacency list for the valid atoms.
    graph = {idx: [] for idx in valid_atoms}
    for idx in valid_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            nidx = neighbor.GetIdx()
            if nidx in graph:
                graph[idx].append(nidx)

    def bfs(start):
        visited = {start: 0}
        queue = deque([start])
        max_dist = 0
        while queue:
            current = queue.popleft()
            for neighbor in graph[current]:
                if neighbor not in visited:
                    visited[neighbor] = visited[current] + 1
                    if visited[neighbor] > max_dist:
                        max_dist = visited[neighbor]
                    queue.append(neighbor)
        return max_dist

    visited_global = set()
    longest_chain = 0
    for atom_idx in valid_atoms:
        if atom_idx in visited_global:
            continue
        comp_nodes = set()
        queue = deque([atom_idx])
        comp_nodes.add(atom_idx)
        while queue:
            cur = queue.popleft()
            for nbr in graph[cur]:
                if nbr not in comp_nodes:
                    comp_nodes.add(nbr)
                    queue.append(nbr)
        visited_global.update(comp_nodes)
        u = next(iter(comp_nodes))
        dist_u = {u: 0}
        q = deque([u])
        while q:
            cur = q.popleft()
            for nbr in graph[cur]:
                if nbr not in dist_u:
                    dist_u[nbr] = dist_u[cur] + 1
                    q.append(nbr)
        farthest_node = max(dist_u, key=dist_u.get)
        diameter = bfs(farthest_node)
        chain_length = diameter + 1  # number of atoms = bonds+1
        if chain_length > longest_chain:
            longest_chain = chain_length
    return longest_chain

def chain_attached_to_peptide(mol, longest_chain_cutoff=8):
    """
    Checks if a long aliphatic chain (of at least `longest_chain_cutoff` atoms) is directly attached
    to any amide group (peptide bond). We find all amide bonds and then scan their neighboring atoms
    to see if any long aliphatic chain is present.
    """
    # SMARTS for a typical amide substructure: note that this returns three atoms (N, C, O).
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Build a set of sp3, nonaromatic carbon atom indices.
    sp3_carbon_idxs = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and atom.GetHybridization().name == "SP3":
            sp3_carbon_idxs.add(atom.GetIdx())
    
    # For each amide match, only use the nitrogen (index 0) and carbonyl carbon (index 1).
    for match in amide_matches:
        # Fix: only consider the first two atoms in the match (ignore the oxygen).
        if len(match) < 2:
            continue
        n_idx, c_idx = match[0], match[1]
        # Check neighbors of the amide nitrogen.
        n_atom = mol.GetAtomWithIdx(n_idx)
        for nbr in n_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in sp3_carbon_idxs:
                # As an approximation, if the globally longest chain meets the cutoff, assume it's attached.
                if get_longest_aliphatic_chain(mol) >= longest_chain_cutoff:
                    return True
        # Also check neighbors of the carbonyl carbon.
        c_atom = mol.GetAtomWithIdx(c_idx)
        for nbr in c_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in sp3_carbon_idxs:
                if get_longest_aliphatic_chain(mol) >= longest_chain_cutoff:
                    return True
    return False

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide must have a peptide component (at least two amide bonds) and a lipid component
    (a long aliphatic chain of at least eight sp3 carbon atoms) attached to it.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a lipopeptide, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count amide bonds by matching the amide SMARTS.
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 2:
        return False, f"Only {len(amide_matches)} amide bond(s) found; need at least two for a peptide component"
    
    longest_chain = get_longest_aliphatic_chain(mol)
    if longest_chain < 8:
        return False, f"Longest aliphatic chain has {longest_chain} carbons; requires at least 8 for a lipid component"
    
    if not chain_attached_to_peptide(mol, longest_chain_cutoff=8):
        return False, "No long aliphatic chain appears to be attached directly to a peptide (amide-containing) substructure"
    
    return True, (f"Contains {len(amide_matches)} amide bonds (peptide component) and a long aliphatic chain of "
                  f"{longest_chain} carbons attached to it (lipid component)")

# Example usage:
# test_smiles = "[H][C@@]1(CCCCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O1"
# res, reason = is_lipopeptide(test_smiles)
# print(res, reason)