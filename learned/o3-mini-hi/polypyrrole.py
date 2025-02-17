"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
Definition: A compound composed of two or more pyrrole units.
A pyrrole unit is defined as a five‐membered aromatic ring containing exactly one nitrogen.
In our improved classifier we not only count pyrrole matches (via SMARTS) but also
require that (a) at least two pyrrole units appear in the same connected substructure and
(b) that the atoms of that cluster form at least 15% of the heavy atoms in the molecule.
"""

from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    
    A polypyrrole is defined as a compound composed of two or more pyrrole units.
    Here a pyrrole unit is detected by a SMARTS pattern for a five-membered aromatic ring
    containing one nitrogen and four carbons. In addition,
    we require that at least two such pyrrole rings are connected (either directly or via adjacent rings)
    and that the combined pyrrole cluster represents a reasonable fraction of the heavy atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polypyrrole, False otherwise.
        str: Reason for the classification.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that matches an aromatic five-membered ring with exactly one nitrogen.
    # [nR] matches an aromatic nitrogen atom in a ring; [cR] matches aromatic carbons in a ring.
    pyrrole_smarts = Chem.MolFromSmarts("[nR]1[cR][cR][cR][cR]1")
    if pyrrole_smarts is None:
        return False, "Error in SMARTS pattern"
    
    # Find all substructure matches (each a tuple of atom indices) matching the pyrrole pattern.
    matches = mol.GetSubstructMatches(pyrrole_smarts, uniquify=True)
    num_pyrroles = len(matches)
    
    if num_pyrroles < 2:
        return False, f"Contains {num_pyrroles} pyrrole unit(s); need at least two for a polypyrrole"
    
    # Deduplicate the matches (each match is a tuple of atom indices)
    unique_matches = [set(match) for match in matches]
    
    # Build connectivity among pyrrole units: two pyrrole matches are 'connected'
    # if any atom in one is directly bonded to any atom in the other.
    n = len(unique_matches)
    # Create a graph as an adjacency list:
    neighbors = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            # Check if any atom in set i is directly bonded to any atom in set j
            connected = False
            for a in unique_matches[i]:
                for b in unique_matches[j]:
                    if mol.GetBondBetweenAtoms(a, b) is not None:
                        connected = True
                        break
                if connected:
                    break
            if connected:
                neighbors[i].add(j)
                neighbors[j].add(i)
    
    # Find connected clusters (using a simple depth-first search)
    visited = set()
    clusters = []
    for i in range(n):
        if i in visited:
            continue
        stack = [i]
        cluster = set()
        while stack:
            current = stack.pop()
            if current not in visited:
                visited.add(current)
                cluster.add(current)
                stack.extend(list(neighbors[current] - visited))
        clusters.append(cluster)
    
    # For each cluster, take the union of atom indices of all pyrrole matches in that cluster.
    cluster_atoms = []
    for cluster in clusters:
        atom_set = set()
        for idx in cluster:
            atom_set |= unique_matches[idx]
        cluster_atoms.append(atom_set)
    
    # Select the largest connected cluster.
    largest_cluster = max(cluster_atoms, key=lambda s: len(s))
    largest_cluster_size = len(largest_cluster)
    
    # Compute total number of heavy atoms (atomic number > 1)
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    
    # Compute the ratio of atoms in the largest pyrrole cluster to all heavy atoms.
    fraction = largest_cluster_size / total_heavy if total_heavy > 0 else 0
    
    # We choose a threshold (e.g., 15%) to require that a significant portion of the molecule comes from connected pyrrole units.
    threshold = 0.15
    if fraction < threshold:
        return False, (f"Found {num_pyrroles} pyrrole unit(s) but the largest connected pyrrole cluster covers "
                       f"{largest_cluster_size} of {total_heavy} heavy atoms ({fraction:.2%}); "
                       "pyrrole units are not sufficiently connected")
    
    return True, (f"Contains {num_pyrroles} pyrrole units with a largest connected cluster covering "
                  f"{largest_cluster_size} of {total_heavy} heavy atoms ({fraction:.2%} of heavy atoms)")

# Example test cases
if __name__ == "__main__":
    test_smiles_list = [
        "c1cc[nH]c1",                           # one pyrrole unit -> not polypyrrole
        "c1cc[nH]c1-c2cc[nH]c2",                  # two pyrrole units connected -> polypyrrole
        "c1cc[nH]c1.Cc2cc[nH]c2",                 # two pyrrole units but unconnected (separate fragments) -> not polypyrrole
    ]
    
    for s in test_smiles_list:
        result, reason = is_polypyrrole(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")