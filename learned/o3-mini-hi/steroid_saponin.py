"""
Classifies: CHEBI:61655 steroid saponin
"""
#!/usr/bin/env python
"""
Classifies: Steroid saponin
Definition: Any saponin derived from a hydroxysteroid.
A steroid saponin is here defined as a molecule that:
 (a) contains a hydroxysteroid nucleus – a fused tetracyclic core (typically three six‐membered rings and one five‐membered ring)
     with a total of roughly 15–20 atoms and at least one hydroxyl (–OH) substituent on that nucleus,
 (b) and contains at least one sugar (glycoside) ring that is connected via an oxygen linker to the nucleus.
Note: Both the nucleus detection and sugar detection use heuristics and may mis‐classify some molecules.
"""
from rdkit import Chem
from itertools import combinations

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): A tuple with a True/False classification and a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Work with explicit hydrogens (to better catch hydroxyls)
    mol = Chem.AddHs(mol)
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    # STEP 1. Identify candidate rings for a steroid nucleus.
    # We require only 5- or 6-membered rings.
    # For a 5-membered candidate, at least 4 atoms should be carbon.
    # For a 6-membered candidate, at least 5 atoms should be carbon.
    candidate_rings = []
    for ring in all_rings:
        if len(ring) not in (5,6):
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if (len(ring) == 5 and carbon_count >= 4) or (len(ring) == 6 and carbon_count >= 5):
            candidate_rings.append((set(ring), len(ring)))
    
    # We need at least four candidate rings to form a tetracyclic (steroid) nucleus.
    if len(candidate_rings) < 4:
        return False, "Not enough candidate rings to form a steroid nucleus"
    
    # STEP 2. Look for 4 candidate rings that likely form a fused steroid nucleus.
    # We expect a nucleus to consist of exactly one 5-membered ring and three 6-membered rings.
    # We relax the union-of-atoms size requirement to between 15 and 20, and require at least 3 pairs of rings share two or more atoms.
    fused_nucleus = None
    for rings4 in combinations(candidate_rings, 4):
        five_count = sum(1 for (_, size) in rings4 if size == 5)
        six_count = sum(1 for (_, size) in rings4 if size == 6)
        if five_count != 1 or six_count != 3:
            continue
        union_atoms = set()
        for (ring_set, _) in rings4:
            union_atoms |= ring_set
        if not (15 <= len(union_atoms) <= 20):
            continue
        # Check fusion: require at least 3 pairs to share at least 2 atoms.
        fusion_pairs = 0
        for (ring1, _), (ring2, _) in combinations(rings4, 2):
            if len(ring1.intersection(ring2)) >= 2:
                fusion_pairs += 1
        if fusion_pairs < 3:
            continue
        
        # This set of rings is our candidate nucleus.
        fused_nucleus = union_atoms
        break

    if fused_nucleus is None:
        return False, "No fused tetracyclic system consistent with a steroid nucleus found"

    # STEP 3. Check that at least one atom in the candidate nucleus bears a hydroxyl (-OH) group.
    hydroxyl_found = False
    for idx in fused_nucleus:
        atom = mol.GetAtomWithIdx(idx)
        # Look at neighbors: if an oxygen is attached and that oxygen has at least one hydrogen, count as hydroxyl.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                    hydroxyl_found = True
                    break
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "Steroid nucleus found, but no hydroxyl substituent detected on it"
    
    # STEP 4. Look for a sugar (glycoside) ring.
    # Heuristics: a sugar ring is defined as a 5- or 6-membered ring with at least 2 oxygen atoms.
    # We also require that the ring is not completely (or almost completely) part of the nucleus.
    # In addition, we look for an oxygen atom in the sugar ring that directly bonds to an atom in the nucleus (implying a glycosidic linkage).
    sugar_found = False
    for ring in all_rings:
        if len(ring) not in (5,6):
            continue
        oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxygen_count < 2:
            continue
        
        # Allow for some overlap with the nucleus, but not complete: we require that NO MORE than 60% of ring atoms come from the nucleus.
        overlap = len(set(ring).intersection(fused_nucleus))
        if overlap > (0.6 * len(ring)):
            continue
        
        # Check for a connecting glycosidic bond:
        linkage_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 8:
                continue  # look only at oxygen atoms in the sugar ring
            # If any neighbor not in this ring is in the candidate nucleus, we accept it as a glycosidic linkage.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in fused_nucleus and nbr.GetIdx() not in ring:
                    linkage_found = True
                    break
            if linkage_found:
                break
        
        if linkage_found:
            sugar_found = True
            break
    
    if not sugar_found:
        return False, "Steroid nucleus found, but no sugar-like ring (glycoside) with a proper linkage detected"
    
    return True, "Molecule contains a fused hydroxysteroid nucleus with an attached sugar (steroid saponin)"

# Example usage (for testing):
# test_smiles = "O(C1C(C2(C(C3C(C4(C(=CC3)CC(O)CC4OC5OC(C(OC6OC(C(O)C(O)C6O)CO)C(O)C5O)C)CC2)C1)C)C(C(O)CCC(C)C)C)C7OC(C(O)C(O)C7O)CO"
# result, reason = is_steroid_saponin(test_smiles)
# print(result, reason)