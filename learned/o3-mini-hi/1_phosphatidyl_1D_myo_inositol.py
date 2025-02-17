"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: 1-phosphatidyl-1D-myo-inositol
Definition:
    A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer 
    and the phosphatidyl group is located at its position 1.
    
Heuristics in this implementation (improved):
1. Check that the molecule contains at least two ester bonds (a C(=O)O substructure).
2. Look through all six-membered rings that consist solely of carbon atoms. For each such ring,
   count substituents: expect exactly five free hydroxyl groups (oxygen with at least one hydrogen)
   and one oxygen that is NOT free – but instead is linked to a phosphorus atom.
3. If a candidate inositol ring is found, follow the oxygen attached to P and then verify that the P atom
   is connected to exactly two ester bonds (i.e. two oxygens that are attached to carbonyl carbons).
This further check helps to ensure that the phosphatidyl group is diacyl and attached at the proper site.
Note: This is a heuristic method that uses substructure and connectivity analyses,
and it does NOT completely guarantee proper stereochemical assignment.
"""

from rdkit import Chem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as 1-phosphatidyl-1D-myo-inositol, else False.
        str: Explanation for the classification result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- Criterion 1: Global check for diacylglycerol ester bonds ---
    # Count all ester bonds (we look for the C(=O)O substructure)
    ester_smarts = "C(=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester bond(s); at least 2 are required for the phosphatidyl moiety"

    # --- Criterion 2: Identify a candidate six-membered inositol ring ---
    # Look for a six-membered ring composed entirely of carbons.
    # In a correct 1-phosphatidyl-1D-myo-inositol headgroup,
    # 5 of the ring carbons should have a free hydroxyl (-OH) substituent,
    # and 1 carbon (the 1-position) should have an oxygen that is not free,
    # but rather is linked to a phosphorus atom (which should be further part of the diacyl chain).
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_inositol = None
    candidate_phospho_O_idx = None  # index of the oxygen substituent on the ring that is linked to P
    for ring in ring_info:
        if len(ring) != 6:
            continue  # only interested in six-membered rings
        # Check that all atoms in the ring are carbon.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue

        free_oh_count = 0
        phospho_sub_count = 0
        phospho_O_candidate = None
        # Loop through each ring carbon.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            neighbors = atom.GetNeighbors()
            # Only consider neighbors that are not in the ring.
            for nbr in neighbors:
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen is bound to a phosphorus atom.
                    attached_to_P = any(nbr2.GetAtomicNum() == 15 for nbr2 in nbr.GetNeighbors() if nbr2.GetIdx() != atom.GetIdx())
                    # Check if oxygen has at least one hydrogen (free hydroxyl).
                    has_H = nbr.GetTotalNumHs() > 0
                    if attached_to_P:
                        phospho_sub_count += 1
                        phospho_O_candidate = nbr
                    elif has_H:
                        free_oh_count += 1
        # For a genuine 1-phosphatidyl-1D-myo-inositol,
        # we require exactly 1 substitution by a phosphate and 5 free -OH groups.
        if phospho_sub_count == 1 and free_oh_count == 5:
            candidate_inositol = ring
            candidate_phospho_O_idx = phospho_O_candidate.GetIdx()
            break

    if candidate_inositol is None:
        return False, "No six-membered inositol ring with 1 phosphate substituent and 5 free hydroxyls found"

    # --- Criterion 3: Check connectivity of the phosphate group ---
    # From the candidate oxygen on the ring (attached to P), follow to the phosphorus atom,
    # and ensure that the P is connected to exactly two ester bonds (i.e. two oxygens that connect to a carbonyl carbon).
    phospho_O = mol.GetAtomWithIdx(candidate_phospho_O_idx)
    phosphorus = None
    for nbr in phospho_O.GetNeighbors():
        if nbr.GetAtomicNum() == 15:
            phosphorus = nbr
            break
    if phosphorus is None:
        return False, "Phospho substituent not attached to a phosphorus atom"

    # Count the number of ester bonds on the phosphorus (excluding the oxygen coming from the ring).
    ester_bond_count = 0
    for nbr in phosphorus.GetNeighbors():
        # Skip the oxygen from the inositol ring.
        if nbr.GetIdx() == phospho_O.GetIdx():
            continue
        if nbr.GetAtomicNum() != 8:
            continue
        # Look for an ester oxygen: attached to a carbon that is doubly bonded to oxygen.
        for nbr2 in nbr.GetNeighbors():
            # Exclude the phosphorus
            if nbr2.GetIdx() == phosphorus.GetIdx():
                continue
            # Check if nbr2 is carbon and has at least one double bond to oxygen.
            if nbr2.GetAtomicNum() == 6:
                bonds = nbr2.GetBonds()
                for bond in bonds:
                    # Check for double bond to oxygen.
                    if bond.GetBondTypeAsDouble() == 2.0:
                        atom1 = bond.GetBeginAtom()
                        atom2 = bond.GetEndAtom()
                        if (atom1.GetAtomicNum() == 8 or atom2.GetAtomicNum() == 8):
                            ester_bond_count += 1
                            break
                break
    if ester_bond_count != 2:
        return False, f"Phosphorus atom is connected to {ester_bond_count} ester bond(s) instead of 2 expected for a diacyl phosphatidyl moiety"

    # If all criteria are met, the molecule is accepted.
    return True, ("Molecule contains a candidate 1D-myo-inositol headgroup (with 1 phospho substitution at the ring "
                  "and 5 free hydroxyls) and its phosphate is linked to two ester bonds as expected for a diacyl phosphatidyl moiety.")


# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with one of the provided correct examples:
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC"
    result, explanation = is_1_phosphatidyl_1D_myo_inositol(test_smiles)
    print("Classification:", result)
    print("Explanation:", explanation)