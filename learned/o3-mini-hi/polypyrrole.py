"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
Definition: A compound composed of two or more pyrrole units.
A pyrrole unit is defined as a five‐membered aromatic ring containing exactly one nitrogen.
In this classifier we:
  - Use the molecule’s ring information (via GetRingInfo) to identify 5‐membered rings.
  - Accept rings that are aromatic and have exactly one nitrogen atom.
  - Group pyrrole rings that are “connected” (if any atom of one ring is either the same or 
    directly bonded to an atom in another ring).
  - Require that one connected cluster contains at least two rings and that the union of atoms 
    in that cluster covers at least 10% of the heavy atoms (atomic number >1) in the molecule.
If the SMILES cannot be parsed or no valid connected cluster is found, the function returns False.
"""

from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Classifies a molecule as a polypyrrole or not.
    
    A polypyrrole is defined as a compound that contains at least two connected pyrrole units.
    Here a pyrrole unit is detected by:
      - Finding a 5-membered ring (via RDKit ring information),
      - Checking that the ring is aromatic (each atom GetIsAromatic() True), and 
      - That the ring contains exactly one nitrogen atom.
    
    Two pyrrole rings are considered connected if at least one atom in one ring is the same as,
    or is directly bonded to, an atom in the other ring.
    
    Finally, the union of atoms belonging to the largest connected pyrrole cluster must represent
    at least 10% of the total heavy atoms (atomic number > 1) in the molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polypyrrole, otherwise False.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # tuple of tuples (each is a set of indices)
    
    pyrrole_rings = []  # will store each pyrrole ring as a frozenset of atom indices
    for ring in all_rings:
        if len(ring) != 5:
            continue  # only consider 5-membered rings
        # Get the atoms for this ring
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check that every atom in the ring is aromatic.
        if not all(atom.GetIsAromatic() for atom in atoms):
            continue
        # Count the nitrogens
        ncount = sum(1 for atom in atoms if atom.GetAtomicNum() == 7)
        if ncount != 1:
            continue
        # Accept this ring as a pyrrole unit.
        pyrrole_rings.append(frozenset(ring))
    
    num_rings = len(pyrrole_rings)
    if num_rings < 2:
        return False, f"Contains {num_rings} pyrrole unit(s); need at least two for a polypyrrole"
    
    # Build connectivity graph between pyrrole rings.
    # Two rings are considered connected if any atom in one ring is either the same as an atom in 
    # the other (i.e. they share an atom) or if any atom in one is directly bonded to any atom in the other.
    n = len(pyrrole_rings)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            connected = False
            for a in pyrrole_rings[i]:
                for b in pyrrole_rings[j]:
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
    
    # Find connected clusters among the pyrrole rings using a depth-first search.
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
    
    # For each cluster, take the union of all atom indices in the pyrrole rings included.
    cluster_atoms_list = []
    for cluster in clusters:
        atoms_in_cluster = set()
        for ring_idx in cluster:
            atoms_in_cluster |= pyrrole_rings[ring_idx]
        cluster_atoms_list.append(atoms_in_cluster)
    
    # Select the largest connected cluster (by atom count)
    largest_cluster_atoms = max(cluster_atoms_list, key=lambda s: len(s))
    lc_size = len(largest_cluster_atoms)

    # Count total heavy atoms in the molecule (atomic number > 1)
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    fraction = lc_size / total_heavy if total_heavy else 0

    # Require that the largest cluster covers at least 10% of the heavy atoms.
    threshold = 0.10
    if fraction < threshold:
        return False, (f"Found {num_rings} pyrrole unit(s) but the largest connected pyrrole cluster covers "
                       f"{lc_size} of {total_heavy} heavy atoms ({fraction:.2%}); "
                       "pyrrole units are not sufficiently connected")
    
    return True, (f"Contains {num_rings} pyrrole units with a largest connected cluster covering "
                  f"{lc_size} of {total_heavy} heavy atoms ({fraction:.2%} of heavy atoms)")

# Example test cases
if __name__ == "__main__":
    test_smiles = [
        "c1cc[nH]c1",                         # one pyrrole unit -> not polypyrrole
        "c1cc[nH]c1-c2cc[nH]c2",                # two pyrrole units connected by a bond -> polypyrrole
        "c1cc[nH]c1.Cc2cc[nH]c2"                # two pyrrole units in separate fragments -> not polypyrrole
    ]
    for s in test_smiles:
        result, reason = is_polypyrrole(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")