"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
Definition: An oxo steroid is defined (heuristically) as a molecule with a fused steroid nucleus
that resembles a cyclopentanoperhydrophenanthrene system (typically 3 six‐membered rings and 1 five‐membered ring,
with roughly 17 carbon atoms in the nucleus) and that bears at least one ketone (C=O) group on a six‐membered ring.
This implementation:
  - Searches for candidate rings (only rings of size 5 or 6).
  - Builds a graph of fused rings (rings sharing at least 2 atoms) and selects the largest connected component.
  - Accepts only if the fused system has 4–5 rings and the carbon count in that nucleus is in the range 15–21.
  - Checks that the nucleus contains at least one five–membered ring and enough six–membered rings (3 if typical, 2 if “nor–steroid”).
  - Finally, looks for a ketone group (C=O) on any six–membered ring in the nucleus.
If any condition fails the molecule is rejected.
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is defined heuristically as a molecule having a fused steroid nucleus
    (a set of 4–5 rings with ~15-21 carbons including at least one five–membered ring)
    and having at least one ketone (C=O) substituent located on a six–membered ring of that nucleus.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a 3-oxo steroid, False otherwise.
       str : Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings (as tuples of atom indices) from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found; not a steroid"

    # Consider only rings of size 5 or 6.
    candidate_rings = [ring for ring in ring_info if len(ring) in (5, 6)]
    if len(candidate_rings) < 4:
        return False, "Less than 4 candidate rings (size 5 or 6); not a typical steroid nucleus"
    
    # Build a graph where each node is a candidate ring.
    # Two rings are connected if they share at least 2 atoms (i.e. are fused).
    n = len(candidate_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if len(set(candidate_rings[i]).intersection(candidate_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components among candidate rings.
    visited = set()
    components = []
    for i in range(n):
        if i not in visited:
            comp = set()
            stack = [i]
            while stack:
                cur = stack.pop()
                if cur not in comp:
                    comp.add(cur)
                    visited.add(cur)
                    stack.extend(ring_graph[cur] - comp)
            components.append(comp)
    
    # Choose the largest connected component as candidate steroid nucleus.
    largest_comp = max(components, key=lambda c: len(c))
    
    # Require that the fused ring system is typical – it must have 4 or 5 rings.
    if len(largest_comp) < 4 or len(largest_comp) > 5:
        return False, "Fused ring system does not have the expected 4 (or 4-5) rings for a steroid nucleus"
    
    # Collect all atoms in the nucleus and count six–membered and five–membered rings.
    nucleus_atoms = set()
    six_count = 0
    five_count = 0
    for idx in largest_comp:
        ring = candidate_rings[idx]
        nucleus_atoms.update(ring)
        if len(ring) == 6:
            six_count += 1
        elif len(ring) == 5:
            five_count += 1

    # Count carbon atoms in the fused nucleus.
    nucleus_carbon_count = sum(1 for atom in mol.GetAtoms() 
                               if atom.GetAtomicNum() == 6 and atom.GetIdx() in nucleus_atoms)
    
    # Check if the carbon count falls in the typical range for steroid nuclei (~17 carbons)
    if nucleus_carbon_count < 15 or nucleus_carbon_count > 21:
        return False, f"Fused ring nucleus has {nucleus_carbon_count} carbons, outside the typical range (15-21)"
    
    # The classic steroid nucleus contains at least one five-membered ring (D-ring)
    if five_count < 1:
        return False, "Fused ring nucleus does not contain a five-membered ring typical of steroids"
    
    # Expect 3 six-membered rings for typical steroids. In nor–steroids (with fewer carbons)
    # we relax this requirement to 2.
    required_six = 3 if nucleus_carbon_count >= 17 else 2
    if six_count < required_six:
        return False, f"Nucleus has only {six_count} six-membered ring(s); expected at least {required_six}"
    
    # Now, search for a ketone (C=O) group on any six-membered ring within the nucleus.
    # Loop over candidate rings from the fused nucleus that are six-membered.
    for idx in largest_comp:
        ring = candidate_rings[idx]
        if len(ring) == 6:
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                # Look at carbon atoms as the carbonyl should be on a carbon.
                if atom.GetAtomicNum() == 6:
                    # Check bonds for a double bond to oxygen.
                    for bond in atom.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                # Found a ketone group on a six-membered ring in the nucleus.
                                return True, ("Found ketone group on a six-membered ring of the steroid nucleus "
                                              "(likely at position 3)")
    
    return False, "No ketone group found on any six-membered ring of the fused steroid nucleus"


# Example usage for testing:
if __name__ == "__main__":
    # A few test molecules (SMILES) from the provided examples:
    test_smiles = [
        "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO",  # methyl prednisolone-16alpha-carboxylate (should be True)
        "C[C@]12CC[C@H]3[C@@H](CC=C4CC(=O)CC[C@]34C)[C@@H]1CCC2=O",  # androst-5-ene-3,17-dione (should be True)
        "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC4=CC(=O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCCC(C)C",   # example steroid nucleus variant (should be True)
        "C[C@H](CCC(=O)NCCCC(O)=O)[C@H]1CC[C@H]2[C@H]3[C@H](CC[C@]12C)[C@@]4(C)CC[C@@H](O)C[C@H]4CC3=O"  # false positive candidate (should be False)
    ]
    for sm in test_smiles:
        result, reason = is_3_oxo_steroid(sm)
        print(sm[:70] + " ...", "->", result, "|", reason)