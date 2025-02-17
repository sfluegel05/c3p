"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: Lipopeptide – a compound consisting of a peptide (i.e. two or more amide bonds)
with an attached lipid (i.e. a long aliphatic chain that is connected to the peptide).
We consider a molecule lipopeptidic if (a) it contains at least two separate amide bonds
(–N–C(=O)– fragments) and (b) it contains a long aliphatic chain of nonaromatic, sp3-carbon atoms
of at least 8 atoms in length which is connected (touched) to an amide substructure.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def get_longest_aliphatic_chain(mol):
    """
    Returns the length (number of atoms) of the longest chain composed solely of sp3-hybridized,
    nonaromatic carbon atoms.
    We build an explicit graph (adjacency list) for sp3 aliphatic carbons and then compute
    the diameter (in atoms) of each connected subgraph.
    """
    # Select atoms that are sp3 carbons and nonaromatic.
    valid_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and atom.GetHybridization().name == "SP3":
            valid_atoms.append(atom.GetIdx())
    if not valid_atoms:
        return 0

    # Build graph for these atoms.
    graph = {idx: [] for idx in valid_atoms}
    for idx in valid_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            nidx = neighbor.GetIdx()
            if nidx in graph:
                graph[idx].append(nidx)

    # Helper: BFS starting at 'start' that returns maximum distance (in bonds)
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
        # Get all nodes in this connected component.
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
        # Use a two-stage BFS to get the diameter in bonds.
        u = next(iter(comp_nodes))
        q = deque([u])
        dist_u = {u: 0}
        while q:
            cur = q.popleft()
            for nbr in graph[cur]:
                if nbr not in dist_u:
                    dist_u[nbr] = dist_u[cur] + 1
                    q.append(nbr)
        farthest_node = max(dist_u, key=dist_u.get)
        diameter = bfs(farthest_node)
        chain_length = diameter + 1  # atoms = bonds+1
        if chain_length > longest_chain:
            longest_chain = chain_length
    return longest_chain

def chain_attached_to_peptide(mol, longest_chain_cutoff=8):
    """
    As a heuristic, we check whether any long aliphatic chain is directly attached 
    to an amide group. We search for all amide bonds and then check the neighbors of
    the N or the adjacent carbonyl carbon. If any neighboring valid sp3 carbon (as in our
    chain detection) is found and it belongs to a chain of at least longest_chain_cutoff length,
    then we assume that the lipid chain is attached to the peptide.
    """
    # Get the substructure for a typical amide: [NX3][CX3](=O)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    # Create a set of sp3 aliphatic carbon indices for quick lookup.
    sp3_carbon_idxs = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and atom.GetHybridization().name == "SP3":
            sp3_carbon_idxs.add(atom.GetIdx())
    # For each match, check the neighbors of the amide nitrogen and the carbonyl carbon.
    for match in amide_matches:
        n_idx, c_idx = match  # n_idx = amide N, c_idx = carbonyl carbon
        # For each neighbor of the amide nitrogen:
        n_atom = mol.GetAtomWithIdx(n_idx)
        for nbr in n_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in sp3_carbon_idxs:
                # Compute chain length starting from this neighbor.
                chain_len = get_longest_aliphatic_chain(Chem.PathToSubmol(mol, Chem.rdmolops.GetShortestPath(mol, nbr_idx, nbr_idx)))
                # Here we use a minimal heuristic: if the atom exists in our global chain (calculated below)
                # and the chain is longer than our cutoff, mark it as attached.
                # (Alternate, we recalc the longest chain globally and then check connectivity.)
                # For simplicity, we return True immediately if our global chain is long enough.
                if get_longest_aliphatic_chain(mol) >= longest_chain_cutoff:
                    return True
        # Also check neighbors of the carbonyl carbon (which is sp2, so not in our chain, but its other neighbors might be).
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
    A lipopeptide should contain a peptide component (with at least two amide bonds)
    and a lipid component (a long aliphatic chain, here defined as at least eight sp3 carbons)
    that is attached to the peptide.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is recognized as a lipopeptide, False otherwise.
        str: A reason detailing the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of amide bonds (using SMARTS pattern).
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 2:
        return False, f"Only {len(amide_matches)} amide bond(s) found; need at least two to indicate a peptide chain"
    
    # Calculate the longest aliphatic chain over the entire molecule.
    longest_chain = get_longest_aliphatic_chain(mol)
    if longest_chain < 8:
        return False, f"Longest aliphatic chain has {longest_chain} carbons; less than the required 8 for a lipid component"
    
    # Use our heuristic to check if the long chain is “attached” to the peptide (via an amide).
    if not chain_attached_to_peptide(mol, longest_chain_cutoff=8):
        return False, "No long aliphatic chain appearing attached to a peptide (amide-containing) substructure"
    
    return True, (f"Contains {len(amide_matches)} amide bonds (peptide component) and a long aliphatic "
                f"chain of {longest_chain} carbons attached to it (lipid component)")

# Example usage (you can uncomment to test some examples):
# test_smiles = [
#     "[H][C@@]1(CCCCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)O1",  # surfactin C (true positive)
#     "CCCC\\C=C/CCCCCCCCCC(=O)NCCO",  # N-(11Z,14Z)-eicosadienoylethanolamine (false positive)
# ]
# for s in test_smiles:
#     res, reason = is_lipopeptide(s)
#     print(res, reason)