"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
Definition: Any oxo steroid where a ketone (oxo) substituent is located at position 3.
Heuristic:
  1. The molecule must have a fused tetracyclic ring system with exactly 3 six‐membered rings 
     and 1 five‐membered ring (the typical natural steroid nucleus).
  2. At least one ketone group (a carbon with a double bond to an oxygen) must be located 
     on one of the six‐membered rings of this fused system.
If no match is found or if the fused system does not match a typical steroid, then the molecule is rejected.
"""

from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is defined (heuristically) as a molecule with a steroid nucleus –
    a fused tetracyclic ring system with (typically) 3 six-membered rings and 1 five-membered ring –
    AND having a ketone (C=O) substituent located on one of the six-membered rings (likely position 3).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3-oxo steroid, False otherwise.
        str : A reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (each ring as a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found, not a steroid"
    
    # Consider rings of size 5 or 6 as candidates
    candidate_rings = [ring for ring in ring_info if len(ring) in (5, 6)]
    if len(candidate_rings) < 4:
        return False, "Less than 4 rings (size 5 or 6) found, not a typical steroid nucleus"
    
    # Build a graph of fused rings:
    # Two candidate rings are considered fused if they share 2 or more atoms.
    n = len(candidate_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(set(candidate_rings[i]).intersection(candidate_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find the largest connected component of fused candidate rings.
    visited = set()
    components = []
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                current = stack.pop()
                if current not in comp:
                    comp.add(current)
                    visited.add(current)
                    stack.extend(ring_graph[current] - comp)
            components.append(comp)
    
    largest_comp = max(components, key=lambda c: len(c))
    if len(largest_comp) < 4:
        return False, "No fused ring system of at least 4 rings detected, not a steroid nucleus"
    
    # Now narrow our attention to the rings in this largest fused system.
    # Count how many six-membered and five-membered rings are present.
    six_count = 0
    five_count = 0
    six_membered_atoms = set()
    for idx in largest_comp:
        ring = candidate_rings[idx]
        if len(ring) == 6:
            six_count += 1
            six_membered_atoms.update(ring)
        elif len(ring) == 5:
            five_count += 1
    # For a typical steroid nucleus we require exactly 3 six-membered rings and 1 five-membered ring.
    if six_count != 3 or five_count != 1:
        return False, f"Fused ring system detected but does not match steroid pattern (found {six_count} six-membered and {five_count} five-membered rings)"
    
    # Now, search for a ketone group (C=O) located on one of the six-membered rings of the steroid nucleus.
    # We iterate over all atoms; if an atom is carbon, is in a ring, belongs to our fused six-membered rings,
    # and is double-bonded to an oxygen atom then we consider it our candidate ketone.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() in six_membered_atoms:
            # Check for a double bond to an oxygen.
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8:
                        return True, "Found ketone group on a six-membered ring of the steroid nucleus (likely at position 3)"
    
    return False, "No ketone group on a six-membered ring of the steroid nucleus, not a 3-oxo steroid"

# Example usage (for testing; uncomment these lines to test):
# test_smiles = [
#     "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO",  # methyl prednisolone-16alpha-carboxylate
#     "C[C@]12CC[C@H]3[C@@H](CC=C4CC(=O)CC[C@]34C)[C@@H]1CCC2=O",  # androst-5-ene-3,17-dione
#     "C1=CC(N[C@]2([C@]1([C@@]3([C@@](CC2)([C@]4([C@](CC3)([C@](CC4)(C(NC(C)(C)C)=O)[H])C)[H])[H])[H])C)[H])=O",  # finasteride (false negative previously)
# ]
# for sm in test_smiles:
#     result, reason = is_3_oxo_steroid(sm)
#     print(sm[:50], "->", result, "|", reason)