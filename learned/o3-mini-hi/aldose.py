"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: Aldose – Aldehydic parent sugars (polyhydroxy aldehydes H[CH(OH)]nC(=O)H, n>=2) 
and their intramolecular hemiacetals.

This improved implementation uses several heuristic rules:
1. The molecule must contain only H, C and O and have 4–8 carbon atoms.
2. Molecules having any internal carbonyl group (i.e. a ketone) are rejected.
3. For an open–chain aldose, we require:
     • exactly one aldehyde pattern ([CH1](=O))
     • that aldehyde carbon is terminal (has a single heavy-atom neighbor)
     • the number of free hydroxyl (-OH) groups (detected by [OX2H]) equals (n_carbons - 1)
4. For a cyclic aldose, we require the existence of one 5– or 6–membered ring that:
     • has exactly one ring oxygen,
     • none of its ring carbons is involved in a carbonyl bond,
     • and the overall number of free –OH groups is (n_carbons - 1)
5. If none of the above criteria are met, the molecule is rejected.

Note: This heuristic does not catch every edge case; sugar chemistry is nuanced,
and small changes in the structure may cause the detection to succeed or fail.
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is a parent aldose (open-chain or cyclic hemiacetal) based on its SMILES string.

    The improved heuristic:
      - Only allows C, H and O.
      - Requires 4 to 8 carbon atoms.
      - Rejects molecules that contain carbonyl groups other than a terminal aldehyde.
      - For an open-chain aldose, expects one terminal aldehyde group and (n-1) free hydroxyl groups.
      - For a cyclic aldose, expects a 5- or 6-membered ring with exactly one oxygen (the ring oxygen)
        and no carbonyl on ring carbons, and the free hydroxyl count equals (n-1).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): Tuple of (True, reason) if the molecule is classified as an aldose;
                     (False, reason) if not;
                     (None, None) if the classification cannot be done.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Rule A: Only allow H, C and O.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, f"Molecule contains {atom.GetSymbol()}, not typical for a parent aldose"
    
    # Count the number of carbons.
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (4 <= n_carbons <= 8):
        return False, f"Carbon count ({n_carbons}) is not in the 4-8 range for a parent aldose"
    
    # Define SMARTS patterns.
    aldehyde_pat = Chem.MolFromSmarts("[CH1](=O)")
    oh_pat = Chem.MolFromSmarts("[OX2H]")  # matches hydroxyl groups

    # Check for non-aldehydic carbonyl groups. In a proper parent aldose the only carbonyl allowed is the terminal aldehyde.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # consider carbon atoms
            # Look for double bonds to oxygen:
            for neighbor in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if neighbor.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                    # Determine if this is an aldehyde carbon: it should be CH1.
                    # (Even though hydrogens are implicit, we can check the explicit degree.)
                    # If not, reject.
                    num_explicit_H = atom.GetTotalNumHs()  # includes implicit H count
                    if num_explicit_H != 1:
                        return False, "Non-terminal carbonyl (ketone) group detected, not an aldose"
    
    # Count free hydroxyl groups.
    oh_matches = mol.GetSubstructMatches(oh_pat)
    n_oh = len(oh_matches)
    
    # Helper: check if free OH count equals expected for a parent aldose.
    expected_oh = n_carbons - 1  # For a parent aldose this is expected.
    
    # Try to find a cyclic candidate.
    ring_info = mol.GetRingInfo()
    candidate_ring = None
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6):
            continue
        # Count oxygens in the ring and record ring carbons.
        ring_oxygen_count = 0
        ring_carbon_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                ring_carbon_indices.append(idx)
        # A typical pyranose/furanose ring has exactly one oxygen.
        if ring_oxygen_count != 1:
            continue
        # Ensure none of the ring carbons is carbonyl-bound.
        found_carbonyl = False
        for idx in ring_carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                    found_carbonyl = True
                    break
            if found_carbonyl:
                break
        if found_carbonyl:
            continue
        candidate_ring = ring
        break

    if candidate_ring is not None:
        # For cyclic aldoses we expect the free OH count to be (n_carbons - 1).
        if n_oh == expected_oh:
            return True, "Cyclic hemiacetal structure consistent with an aldose detected (sugar ring with proper hydroxyl pattern)"
        else:
            return False, f"Cyclic candidate ring found but hydroxyl count ({n_oh}) does not match expected ({expected_oh}) for a parent aldose"
    
    # Otherwise try open-chain classification.
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pat)
    if len(aldehyde_matches) != 1:
        return False, f"Expected exactly one terminal aldehyde group in open-chain aldose, found {len(aldehyde_matches)}"
    # Check that the aldehyde carbon is terminal (only one heavy-atom neighbor).
    aldehyde_idx = aldehyde_matches[0][0]
    aldehyde_atom = mol.GetAtomWithIdx(aldehyde_idx)
    heavy_neighbors = [nbr for nbr in aldehyde_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
    if len(heavy_neighbors) != 1:
        return False, "Aldehyde group is not terminal"
    if n_oh != expected_oh:
        return False, f"Open-chain aldehyde detected but hydroxyl count ({n_oh}) does not match expected ({expected_oh})"
    return True, "Open-chain molecule with terminal aldehyde and correct hydroxyl pattern detected"

# Example test calls:
if __name__ == '__main__':
    test_smiles = [
        # True positives (cyclic aldoses and open-chain aldoses):
        "O[C@H]1COC(O)[C@H](O)[C@@H]1O",             # L-lyxopyranose (cyclic)
        "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",      # D-allopyranose (cyclic)
        "OC[C@@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O", # alpha-L-gulose (cyclic)
        "OC[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1O",    # alpha-D-gulose (cyclic)
        "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O",    # beta-D-gulose (cyclic)
        "[C@@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)CO",# beta-L-idopyranose (cyclic)
        "[H]C([H])([C@]([H])(O)C=O)[C@]([H])(O)[C@@]([H])(C)O", # tyvelose (open-chain candidate)
        "OC[C@@H](O)[C@@H](O)C=O",                      # D-erythrose (open-chain candidate)
        "O1[C@H](O)[C@H](O)[C@H](O)[C@@H](O)C1",         # beta-L-lyxopyranose (cyclic)
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",     # alpha-D-glucose (cyclic)
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",   # aldehydo-D-galactose (open-chain candidate)
        # The following examples were previously false positives:
        "OC[C@H]1O[C@](O)(CO)[C@@H](O)[C@H](O)[C@@H]1O", # beta-D-sedoheptulopyranose (should be ketose so rejected)
        "O1C(OCCO)C(O)C(O)C(O)C1C(O)=O",                # 3,4,5-trihydroxy...oxane-2-carboxylic acid (contains extra carbonyl)
        "OC[C@H]1OC(O)(CO)[C@@H](O)[C@H]1O",             # D-tagatofuranose (ketose, not aldose)
        # Previously false negatives:
        "[H]C([H])([C@]([H])(O)C=O)[C@]([H])(O)[C@@]([H])(C)O", # tyvelose (should now be accepted if free OH count fits)
        "OC[C@@H](O)[C@@H](O)C=O",                       # D-erythrose
        "O1[C@H]([C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O)C",   # 6-deoxy-beta-L-talopyranose (deoxy: free OH count off)
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",  # aldehydo-D-galactose (if aldehyde not terminal, will be rejected)
        "[H]C(=O)[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO", # aldehydo-L-gulose (aldehyde not terminal)
        "[H][C@](C)(O)[C@@]([H])(O)[C@]([H])(O)[C@]([H])(O)C=O", # aldehydo-D-rhamnose (deoxy pattern)
        "O1[C@@H]([C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O)C", # 6-deoxy-beta-D-gulopyranose (deoxy: free OH count off)
        "C[C@@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O",           # L-rhamnopyranose (deoxy)
        "C[C@H]1O[C@H](O)[C@H](O)C[C@H]1O",              # alpha-abequopyranose (deoxy pattern)
        "O1[C@H]([C@H](O)[C@@H](O)C1O)[C@H](O)C",         # D-fucofuranose (ketose pattern)
    ]
    for s in test_smiles:
        result, reason = is_aldose(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*50}")