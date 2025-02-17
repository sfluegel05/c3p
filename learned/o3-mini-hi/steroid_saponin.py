"""
Classifies: CHEBI:61655 steroid saponin
"""
#!/usr/bin/env python
"""
Classifies: Steroid saponin
Definition: Any saponin derived from a hydroxysteroid.
A steroid saponin is defined here as a molecule that (a) contains a hydroxysteroid nucleus—that is,
a fused tetracyclic core (approximately 15–19 atoms in total) composed of 4 rings (typically 3 six‐membered and 1 five‐membered)
in which at least one atom of the core bears a hydroxyl (–OH) substituent—and (b) contains at least one sugar moiety.
In this implementation we try to:
  • Find candidate rings (5‐ or 6‐membered rings with mostly carbon atoms) that could belong to a steroid nucleus.
  • Search among combinations of 4 candidate rings for one that:
      - contains exactly one 5‐membered and three 6‐membered rings,
      - the union of their atoms is between 16 and 19,
      - and they show evidence of “fusion” by sharing two or more atoms in at least 3 distinct ring–ring pairs.
  • Check at least one atom in that nucleus is substituted by a hydroxyl.
  • Look for the presence of at least one sugar ring: defined as a 5‐ or 6‐membered ring with ≥2 oxygen atoms,
    and not completely overlapping the steroid nucleus.
Note: This heuristic does not guarantee 100% sensitivity/specificity.
"""

from rdkit import Chem
from itertools import combinations

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.

    Args:
      smiles (str): SMILES string of the molecule.

    Returns:
      (bool, str): Tuple of classification (True/False) and a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens
    mol = Chem.AddHs(mol)
    
    # Get ring information from molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # tuple of atom index tuples
    
    # Step 1. Identify candidate rings for the steroid nucleus.
    # A candidate ring must be of size 5 or 6.
    # For a 5-membered ring, require at least 4 atoms to be carbons;
    # for 6-membered ring, require at least 5 carbons.
    candidate_rings = []
    for ring in all_rings:
        if len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if (len(ring) == 5 and carbon_count >= 4) or (len(ring) == 6 and carbon_count >= 5):
            # Save as a tuple (set_of_atom_indices, ring_size)
            candidate_rings.append((set(ring), len(ring)))
    
    # If we do not have at least 4 candidate rings, we cannot form the steroid nucleus.
    if len(candidate_rings) < 4:
        return False, "Not enough candidate rings to form steroid nucleus"
    
    # Step 2. Look for a set of 4 rings that meet our conditions
    fused_nucleus = None  # will hold the union of atom indices if found
    for rings4 in combinations(candidate_rings, 4):
        # Count how many 5-membered and 6-membered rings:
        five_count = sum(1 for (_, size) in rings4 if size == 5)
        six_count = sum(1 for (_, size) in rings4 if size == 6)
        if five_count != 1 or six_count != 3:
            continue
        # Get union of atoms in these rings:
        union_atoms = set()
        for (ring_set, _) in rings4:
            union_atoms |= ring_set
        # Heuristic: Steroid nucleus should have between 16 and 19 distinct atoms
        if not (16 <= len(union_atoms) <= 19):
            continue
        # Check fusion: count how many pairs share at least 2 atoms.
        fusion_pairs = 0
        for (ring1, _), (ring2, _) in combinations(rings4, 2):
            if len(ring1.intersection(ring2)) >= 2:
                fusion_pairs += 1
        # In a typical steroid, there are 3 bonds linking the four rings.
        if fusion_pairs < 3:
            continue
        
        # Found a candidate steroid nucleus.
        fused_nucleus = union_atoms
        break
    
    if fused_nucleus is None:
        return False, "No fused tetracyclic system consistent with a steroid nucleus found"
    
    # Step 3. Check that at least one atom in the candidate nucleus bears a hydroxyl (-OH) substituent.
    hydroxyl_found = False
    for atom_idx in fused_nucleus:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Look through neighbors – if oxygen is attached and that oxygen is in turn bonded to hydrogen.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                    hydroxyl_found = True
                    break
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "Steroid nucleus found, but no hydroxyl substituent detected on it"
    
    # Step 4. Look for a sugar moiety.
    # Define sugar-like ring as a 5- or 6-membered ring (from any of all rings)
    # containing at least 2 oxygen atoms. We also require that this sugar ring is not entirely
    # part of the steroid nucleus (so that it represents an appended glycoside).
    sugar_found = False
    for ring in all_rings:
        if len(ring) not in (5, 6):
            continue
        # Count oxygen atoms in the ring.
        oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxygen_count < 2:
            continue
        # Check overlap with steroid nucleus: if more than half the atoms are part of nucleus,
        # skip (they may be counted as part of the steroid core rather than a glycoside).
        overlap = len(set(ring).intersection(fused_nucleus))
        if overlap > (len(ring) // 2):
            continue
        # Otherwise, we consider it a sugar ring.
        sugar_found = True
        break
    if not sugar_found:
        return False, "Steroid nucleus found, but no sugar-like ring (glycoside) detected"
    
    return True, "Molecule contains a fused hydroxysteroid nucleus with a sugar moiety (steroid saponin)"

# Example usage (for testing):
# test_smiles = "O(C1C(C2(C(C3C(C4(C(=CC3)CC(O)CC4OC5OC(C(OC6OC(C(O)C(O)C6O)CO)C(O)C5O)C)CC2)C1)C)C(C(O)CCC(C)C)C)C7OC(C(O)C(O)C7O)CO"
# result, reason = is_steroid_saponin(test_smiles)
# print(result, reason)