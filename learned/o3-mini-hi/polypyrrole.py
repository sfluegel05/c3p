"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
Definition: A compound composed of two or more pyrrole units.
A pyrrole unit is defined as a five‐membered aromatic ring containing exactly one nitrogen (as –NH)
and four carbons.
This classifier:
  - Uses a SMARTS query to detect pyrrole units.
  - For each match, it verifies that the five ring atoms really are “pure” (i.e. four carbons and one –NH).
  - It then groups pyrrole rings into connected clusters. Two pyrrole units are “connected” if at least one atom from one unit 
    is either the same as or is directly bonded to an atom in the other unit.
  - Finally, the union of atoms in the largest pyrrole cluster must cover at least 10% of the total heavy atoms (atomic number >1).
If the SMILES cannot be parsed or no valid connected pyrrole cluster is found, the function returns (False, reason).
Note: This method is heuristic and may still mis‐classify some edge cases.
"""

from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole.
    
    A polypyrrole must contain at least two connected pyrrole units.
    A pyrrole unit is detected by a SMARTS search for a five-membered aromatic ring
    with exactly four carbons and one nitrogen carrying a hydrogen ([nH]).
    
    Two pyrrole rings are considered connected if any atom in one is identical to,
    or bonded to, an atom in the other.
    
    Finally, the union of atoms present in the largest connected pyrrole cluster must represent 
    at least 10% of all heavy atoms (atomic number > 1) in the molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polypyrrole, otherwise False.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for a pyrrole ring: a 5‐membered aromatic ring with exactly one nitrogen that has an explicit hydrogen.
    # Note: 'nH' matches an aromatic nitrogen with at least one H; the other four positions must be aromatic carbons.
    pyrrole_smarts = "[nH]1cccc1"
    pyrrole_query = Chem.MolFromSmarts(pyrrole_smarts)
    if pyrrole_query is None:
        return False, "Error in pyrrole SMARTS"
    
    # Find pyrrole substructure matches. Each match is a 5-tuple of atom indices.
    matches = mol.GetSubstructMatches(pyrrole_query, uniquify=True)
    
    # Filter matches: require exactly 5 atoms and verify element types.
    pyrrole_units = []
    for match in matches:
        if len(match) != 5:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in match]
        # Count nitrogen atoms and carbon atoms in the match.
        n_atoms = [atom for atom in atoms if atom.GetAtomicNum() == 7]
        c_atoms = [atom for atom in atoms if atom.GetAtomicNum() == 6]
        # A valid pyrrole unit should have exactly 1 nitrogen (and it should be –NH) and 4 carbons.
        # (We assume substituents on the ring are not part of the 5-atom match.)
        if len(n_atoms) != 1 or len(c_atoms) != 4:
            continue
        pyrrole_units.append(frozenset(match))
        
    num_units = len(pyrrole_units)
    if num_units < 2:
        return False, f"Contains {num_units} pyrrole unit(s); need at least two for a polypyrrole"
    
    # Build connectivity graph among the pyrrole unit matches.
    # Two pyrrole units are "connected" if any atom from one is the same as (shared) or
    # is directly bonded (via a single bond) to any atom from the other.
    n = len(pyrrole_units)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
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
    
    # Find connected clusters using depth-first search.
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
    
    # For each cluster, take the union of atom indices from the pyrrole units in that cluster.
    cluster_atoms_list = []
    for cluster in clusters:
        atoms_in_cluster = set()
        for ring_idx in cluster:
            atoms_in_cluster |= pyrrole_units[ring_idx]
        cluster_atoms_list.append(atoms_in_cluster)
    
    # Select the largest cluster (by atom count).
    largest_cluster_atoms = max(cluster_atoms_list, key=lambda s: len(s))
    lc_size = len(largest_cluster_atoms)
    
    # Count heavy atoms (atomic number > 1) in the whole molecule.
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    fraction = lc_size / total_heavy if total_heavy else 0
    
    # Require that the largest connected pyrrole cluster covers at least 10% of heavy atoms.
    threshold = 0.10
    if fraction < threshold:
        return False, (f"Found {num_units} pyrrole unit(s) but the largest connected pyrrole cluster covers "
                       f"{lc_size} of {total_heavy} heavy atoms ({fraction:.2%}); "
                       "pyrrole units are not sufficiently connected")
    
    return True, (f"Contains {num_units} pyrrole unit(s) with a largest connected cluster covering "
                  f"{lc_size} of {total_heavy} heavy atoms ({fraction:.2%} of heavy atoms)")

# Example test cases
if __name__ == "__main__":
    test_smiles = [
        "c1cc[nH]c1",                         # single pyrrole -> not polypyrrole
        "c1cc[nH]c1-c2cc[nH]c2",                # two pyrrole rings directly connected -> polypyrrole
        "c1cc[nH]c1.Cc2cc[nH]c2"                # two pyrrole rings in separate fragments -> not polypyrrole
    ]
    for s in test_smiles:
        result, reason = is_polypyrrole(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")