"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
Definition: A compound composed of two or more pyrrole units.
A pyrrole unit is defined as a five‐membered aromatic ring containing exactly one nitrogen.
This improved classifier not only finds pyrrole rings via SMARTS but also (a) filters only those rings 
whose unique nitrogen bears exactly one hydrogen (as expected for pyrrole) and (b) requires that at least 
two such pyrrole rings are connected (directly or via an adjoining bond) and that the union of atoms 
in that connected cluster comprises at least 10% of the heavy atoms of the molecule.
"""

from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    
    A polypyrrole is defined as a compound composed of two or more pyrrole units.
    Here a pyrrole unit is detected by a SMARTS pattern for a five‐membered aromatic ring
    containing one nitrogen (with exactly one hydrogen). In addition,
    we require that at least two such pyrrole rings are connected (in one cluster)
    and that the combined pyrrole cluster represents at least 10% of the heavy atoms in the molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polypyrrole, False otherwise.
        str: Reason for classification.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that matches a five-membered aromatic ring with one nitrogen.
    pyrrole_smarts = Chem.MolFromSmarts("[nR]1[cR][cR][cR][cR]1")
    if pyrrole_smarts is None:
        return False, "Error in SMARTS pattern"
    
    # Get all ring matches using the SMARTS.
    raw_matches = mol.GetSubstructMatches(pyrrole_smarts, uniquify=True)
    
    # Filter matches to ensure that the single nitrogen in the five-membered ring has exactly one hydrogen.
    filtered_matches = []
    for match in raw_matches:
        atom_indices = list(match)
        # Identify the nitrogen atom(s) in the match
        n_atoms = [a for a in atom_indices if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
        # A proper pyrrole should have exactly one nitrogen.
        if len(n_atoms) != 1:
            continue
        n_atom = mol.GetAtomWithIdx(n_atoms[0])
        # Check that the nitrogen atom has exactly one attached hydrogen (using total H count).
        if n_atom.GetTotalNumHs() != 1:
            continue
        filtered_matches.append(set(match))
    
    num_pyrroles = len(filtered_matches)
    if num_pyrroles < 2:
        return False, f"Contains {num_pyrroles} pyrrole unit(s); need at least two for a polypyrrole"
    
    # Build connectivity between pyrrole matches: two pyrrole rings are connected if any atom in one is directly bonded to an atom in the other.
    n = len(filtered_matches)
    neighbors = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            connected = False
            for a in filtered_matches[i]:
                for b in filtered_matches[j]:
                    if mol.GetBondBetweenAtoms(a, b) is not None:
                        connected = True
                        break
                if connected:
                    break
            if connected:
                neighbors[i].add(j)
                neighbors[j].add(i)
    
    # Find connected clusters among the pyrrole matches using depth-first search.
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
    
    # For each cluster, take the union of all atom indices from the pyrrole matches in that cluster.
    cluster_atoms = []
    for cluster in clusters:
        atoms_in_cluster = set()
        for idx in cluster:
            atoms_in_cluster |= filtered_matches[idx]
        cluster_atoms.append(atoms_in_cluster)
    
    # Select the largest connected cluster.
    largest_cluster = max(cluster_atoms, key=lambda s: len(s))
    largest_cluster_size = len(largest_cluster)
    
    # Count total heavy atoms in the molecule (atomic number > 1)
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    fraction = largest_cluster_size / total_heavy if total_heavy > 0 else 0
    
    # Use a threshold of 10% (lowered from 15%) to allow for true positives that have smaller pyrrole portions.
    threshold = 0.10
    if fraction < threshold:
        return False, (f"Found {num_pyrroles} pyrrole unit(s) but the largest connected pyrrole cluster covers "
                       f"{largest_cluster_size} of {total_heavy} heavy atoms ({fraction:.2%}); "
                       "pyrrole units are not sufficiently connected")
    
    return True, (f"Contains {num_pyrroles} pyrrole units with a largest connected cluster covering "
                  f"{largest_cluster_size} of {total_heavy} heavy atoms ({fraction:.2%} of heavy atoms)")

# Example test cases
if __name__ == "__main__":
    test_smiles = [
        "c1cc[nH]c1",                           # one pyrrole unit -> not polypyrrole
        "c1cc[nH]c1-c2cc[nH]c2",                  # two pyrrole units connected by a bond -> polypyrrole
        "c1cc[nH]c1.Cc2cc[nH]c2"                  # two pyrrole units but in separate fragments -> not polypyrrole
    ]
    for s in test_smiles:
        result, reason = is_polypyrrole(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")