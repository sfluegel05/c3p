"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
Definition: A compound composed of two or more pyrrole units.
A pyrrole unit (for our purposes) is a five‐membered aromatic ring containing exactly four carbons and one nitrogen.
This classifier:
  - Uses RDKit’s ring information to detect aromatic 5‐membered rings.
  - Verifies that exactly one of the ring atoms is nitrogen (and the rest are carbons).
  - Clusters pyrrole rings that are “connected” (i.e. share an atom or are directly bonded).
  - Finally, requires that a cluster of at least two pyrrole rings exists and that 
    the union of atoms in that cluster covers at least 5% of the molecule’s heavy atoms.
Note: This method is heuristic and may mis‐classify some edge cases.
"""

from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole.
    
    A polypyrrole is defined as a compound that contains at least two connected pyrrole units.
    A pyrrole unit is an aromatic five‐membered ring with exactly one nitrogen and four carbons.
    We “cluster” such rings: two pyrrole rings are connected if any atom of one is either identical
    to or directly bonded to an atom in the other.
    
    In addition, we require that the union of atom indices in the largest connected pyrrole cluster 
    represents at least 5% of all heavy atoms (atomic number > 1) in the molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a polypyrrole, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    pyrrole_units = []  # list of sets of atom indices for valid pyrrole rings
    for ring in ring_info:
        # We only consider five-membered rings
        if len(ring) != 5:
            continue
        # Check that every atom in the ring is aromatic
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Count element types
        num_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        num_carbon = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        # We expect exactly one nitrogen and four carbons
        if num_nitrogen == 1 and num_carbon == 4:
            pyrrole_units.append(set(ring))
    
    num_units = len(pyrrole_units)
    if num_units < 2:
        return False, f"Contains {num_units} pyrrole unit(s); need at least two for a polypyrrole"
    
    # Build connectivity graph among the pyrrole units.
    # Two pyrrole rings are connected if any atom in one is the same as,
    # or is directly bonded to an atom in the other.
    n = len(pyrrole_units)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            connected = False
            for a in pyrrole_units[i]:
                for b in pyrrole_units[j]:
                    if a == b:
                        connected = True
                        break
                    if mol.GetBondBetweenAtoms(a, b) is not None:
                        connected = True
                        break
                if connected:
                    break
            if connected:
                graph[i].add(j)
                graph[j].add(i)
    
    # Use depth-first search to group pyrrole units into connected clusters
    visited = set()
    clusters = []
    for i in range(n):
        if i in visited:
            continue
        stack = [i]
        cluster = set()
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            cluster.add(current)
            stack.extend(graph[current] - visited)
        clusters.append(cluster)
    
    # For each cluster, compute the union of atom indices from its pyrrole rings.
    cluster_atoms_list = []
    for cluster in clusters:
        atoms_in_cluster = set()
        for idx in cluster:
            atoms_in_cluster |= pyrrole_units[idx]
        cluster_atoms_list.append(atoms_in_cluster)
    
    # Select the largest cluster (by number of atoms in the union)
    largest_cluster = max(cluster_atoms_list, key=lambda s: len(s))
    lc_atoms = len(largest_cluster)
    
    # Count heavy atoms (atomic number > 1) in the whole molecule.
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    fraction = lc_atoms / total_heavy if total_heavy else 0.0
    
    # Require that the largest connected pyrrole cluster covers at least 5% of heavy atoms.
    threshold = 0.05
    if fraction < threshold:
        return False, (f"Found {num_units} pyrrole unit(s) but the largest connected cluster covers "
                       f"{lc_atoms} of {total_heavy} heavy atoms ({fraction:.2%}); "
                       "pyrrole units are not sufficiently connected")
    
    # Passed all tests: report the number of pyrrole rings and the connected cluster details.
    return True, (f"Contains {num_units} pyrrole unit(s) with a largest connected cluster covering "
                  f"{lc_atoms} of {total_heavy} heavy atoms ({fraction:.2%} of heavy atoms)")

# Example test cases (feel free to add more)
if __name__ == "__main__":
    test_smiles = [
        "c1cc[nH]c1",                          # single pyrrole -> not polypyrrole
        "c1cc[nH]c1-c2cc[nH]c2",                 # two directly connected pyrrole rings -> polypyrrole
        "c1cc[nH]c1.Cc2cc[nH]c2"                 # two pyrrole rings in separate fragments -> not polypyrrole
    ]
    for s in test_smiles:
        result, reason = is_polypyrrole(s)
        print("SMILES:", s)
        print("Result:", result)
        print("Reason:", reason)
        print()