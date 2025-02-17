"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: Lipopeptide – a compound consisting of a peptide (i.e. one or more amide bonds)
with an attached lipid (i.e. a long aliphatic chain).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def get_longest_aliphatic_chain(mol):
    """
    Returns the length of the longest chain (number of atoms) that is:
      – composed solely of nonaromatic, sp3-hybridized carbon atoms.
    Instead of a recursive DFS, we build an explicit graph (dictionary) of valid atoms
    and use iterative BFS to calculate the longest distance (chain length) in each connected
    component.
    """
    # Build a list of atom indices that are sp3 aliphatic carbons.
    valid_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and atom.GetHybridization().name == "SP3":
            valid_atoms.append(atom.GetIdx())

    if not valid_atoms:
        return 0

    # Build an adjacency list (graph) for valid atoms.
    graph = {idx: [] for idx in valid_atoms}
    for idx in valid_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            nidx = neighbor.GetIdx()
            if nidx in graph:  # neighbor is also valid
                graph[idx].append(nidx)

    # Function to perform BFS from a starting atom index.
    def bfs(start):
        visited = {start: 0}  # atom index -> distance (in bonds)
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

    # Instead of restarting BFS from every atom (which can be expensive),
    # we iterate over connected components. For each component we can use a two-stage BFS:
    # first, pick an arbitrary node and find the furthest node; then, run BFS from that furthest node
    # to get the component's diameter.
    visited_global = set()
    longest_chain = 0
    for atom_idx in valid_atoms:
        if atom_idx in visited_global:
            continue
        # Collect all nodes in this connected component via BFS.
        comp_nodes = set()
        queue = deque([atom_idx])
        comp_nodes.add(atom_idx)
        while queue:
            current = queue.popleft()
            for neighbor in graph[current]:
                if neighbor not in comp_nodes:
                    comp_nodes.add(neighbor)
                    queue.append(neighbor)
        visited_global.update(comp_nodes)
        # Pick an arbitrary node 'u' in the component and perform BFS.
        u = next(iter(comp_nodes))
        # Find the furthest node from u.
        dist_u = {}
        q = deque([u])
        dist_u[u] = 0
        while q:
            cur = q.popleft()
            for nbr in graph[cur]:
                if nbr not in dist_u:
                    dist_u[nbr] = dist_u[cur] + 1
                    q.append(nbr)
        # The farthest node from u:
        farthest_node = max(dist_u, key=dist_u.get)
        # Now BFS from farthest_node to determine diameter (longest chain length in bonds).
        diameter = bfs(farthest_node)
        chain_length = diameter + 1  # number of atoms = bonds + 1
        if chain_length > longest_chain:
            longest_chain = chain_length
    return longest_chain

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide should contain a peptide component (evidenced by amide bonds)
    and a lipid component (evidenced by a long aliphatic chain attached).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is recognized as a lipopeptide, False otherwise.
        str: A reason detailing the classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, check for peptide (amide) bonds.
    # Typical peptide/amide bond fragment pattern: -N-C(=O)-
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    peptide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(peptide_matches) < 1:
        return False, "No amide (peptide) bonds found; peptide component missing"
    
    # Next, try to detect the lipid component by calculating the longest 
    # chain of nonaromatic sp3 carbons.
    longest_chain = get_longest_aliphatic_chain(mol)
    if longest_chain < 8:
        return False, f"Longest aliphatic chain has {longest_chain} carbons; lipid component missing"
    
    # If both a peptide signature and a long lipid chain are found, classify as lipopeptide.
    return True, f"Contains amide bonds (peptide component) and a long aliphatic chain of {longest_chain} carbons (lipid component)"

# Example usage (uncomment below lines to test):
# smiles_examples = [
#     "[H][C@@]1(CCCCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)O1",  # surfactin C
#     "OC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCCC"  # long fatty acid (should not be classified as lipopeptide)
# ]
# for s in smiles_examples:
#     result, reason = is_lipopeptide(s)
#     print(result, reason)