"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: 1-phosphatidyl-1D-myo-inositol
Definition:
    A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer 
    and the phosphatidyl group is located at its position 1.

Heuristics in this implementation:
1. Verify that the molecule has at least two ester bonds (C(=O)O) overall.
2. Look for a six-membered ring of carbon atoms that serves as the inositol headgroup.
   In such a ring we expect 5 of its carbons to bear a free hydroxyl (–OH) and one to
   be substituted with an oxygen which is in turn bonded to a phosphorus.
3. From that inositol oxygen, follow to the phosphorus atom and verify that:
   • there are exactly two ester bonds from phosphorus to acyl chains (using the SMARTS pattern "[#15]-O-C(=O)")
   • there is exactly one phosphoryl (P=O) bond (using the SMARTS pattern "[#15]=O")
Note: This is a heuristic method; there will be edge cases that are not captured.
"""

from rdkit import Chem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    
    Args:
        smiles (str): The SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as 1-phosphatidyl-1D-myo-inositol, else False.
        str: An explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ------ Step 1: Global check for at least two ester bonds ------
    # Look for the generic ester substructure "C(=O)O" anywhere.
    generic_ester_smarts = "C(=O)O"
    generic_ester = Chem.MolFromSmarts(generic_ester_smarts)
    if generic_ester is None:
        return False, "Internal error in ester SMARTS"
    ester_matches = mol.GetSubstructMatches(generic_ester)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester bond(s); at least 2 are required for a diacyl phosphatidyl moiety"
    
    # ------ Step 2: Identify candidate inositol ring ------
    # Search for a six-membered ring composed entirely of carbons.
    rings = mol.GetRingInfo().AtomRings()
    candidate_inositol = None
    candidate_phospho_O_idx = None  # index of the oxygen on the ring that attaches to P
    for ring in rings:
        if len(ring) != 6:
            continue  # only consider six-membered rings
        # Verify that all atoms in the ring are carbons.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        
        free_oh_count = 0
        phospho_sub_count = 0
        phospho_O_candidate = None
        # For each ring carbon, examine substituents not in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip atoms in the ring
                # Look for oxygen substituents.
                if nbr.GetAtomicNum() == 8:
                    # If oxygen is bonded to phosphorus then it is a candidate for the phospho substitution.
                    attached_to_P = any(n.GetAtomicNum() == 15 for n in nbr.GetNeighbors() if n.GetIdx() != atom.GetIdx())
                    # Test if oxygen has hydrogen(s) to be a free hydroxyl.
                    has_H = nbr.GetTotalNumHs() > 0
                    if attached_to_P:
                        phospho_sub_count += 1
                        phospho_O_candidate = nbr
                    elif has_H:
                        free_oh_count += 1
        # Expect exactly 1 oxygen leading to phosphate and 5 free hydroxyls.
        if phospho_sub_count == 1 and free_oh_count == 5:
            candidate_inositol = ring
            candidate_phospho_O_idx = phospho_O_candidate.GetIdx()
            break

    if candidate_inositol is None:
        return False, "No six-membered inositol ring with 1 phosphate substituent and 5 free hydroxyls found"
    
    # ------ Step 3: Examine connectivity of the phosphate group ------
    # From the inositol ring oxygen, find the phosphorus atom.
    phospho_O = mol.GetAtomWithIdx(candidate_phospho_O_idx)
    phosphorus = None
    for nbr in phospho_O.GetNeighbors():
        if nbr.GetAtomicNum() == 15:
            phosphorus = nbr
            break
    if phosphorus is None:
        return False, "The oxygen substituent on the inositol ring is not attached to any phosphorus atom"
    
    p_idx = phosphorus.GetIdx()
    
    # Use SMARTS to count ester bonds: pattern "[#15]-O-C(=O)"
    ester_phosphate_smarts = "[#15]-O-C(=O)"
    ester_phosphate_pattern = Chem.MolFromSmarts(ester_phosphate_smarts)
    if ester_phosphate_pattern is None:
        return False, "Internal error in ester-phosphate SMARTS"
    ester_bonds = mol.GetSubstructMatches(ester_phosphate_pattern)
    ester_count = sum(1 for match in ester_bonds if match[0] == p_idx)
    if ester_count != 2:
        return False, (f"Phosphorus atom is connected to {ester_count} ester bond(s) (via O-C(=O) groups) "
                       "instead of 2 expected for a diacyl phosphatidyl moiety")
    
    # Use SMARTS to count phosphoryl bonds: pattern "[#15]=O"
    phosphoryl_smarts = "[#15]=O"
    phosphoryl_pattern = Chem.MolFromSmarts(phosphoryl_smarts)
    if phosphoryl_pattern is None:
        return False, "Internal error in phosphoryl SMARTS"
    phosphoryl_bonds = mol.GetSubstructMatches(phosphoryl_pattern)
    phosphoryl_count = sum(1 for match in phosphoryl_bonds if match[0] == p_idx)
    if phosphoryl_count != 1:
        return False, (f"Phosphorus atom is connected to {phosphoryl_count} phosphoryl (P=O) bond(s); "
                       "expected exactly 1")
    
    return True, ("Molecule contains an inositol headgroup (six-membered carbon ring with 5 free OH and 1 phosphate substitution) "
                  "and its phosphate group is connected to a phosphoryl oxygen and 2 acyl ester oxygens as expected.")

# ----- For testing purposes -----
if __name__ == "__main__":
    # Use one of the provided correct examples for testing.
    test_smiles = ("CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)"
                   "[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC")
    result, explanation = is_1_phosphatidyl_1D_myo_inositol(test_smiles)
    print("Classification:", result)
    print("Explanation:", explanation)