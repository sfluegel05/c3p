"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: 1-phosphatidyl-1D-myo-inositol
Definition:
    A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer 
    and the phosphatidyl group is located at its position 1.
    
Heuristics in this implementation:
1. Verify that the molecule has at least two ester bonds (the C(=O)O substructure).
2. Identify a six-membered inositol ring (a ring of six carbons) that has five free –OH groups and one –OH that is substituted with a phosphorus.
3. From that candidate oxygen substituent on the ring, follow to the phosphorus atom and verify that the phosphate has exactly 4 oxygen neighbors: one from the inositol ring, one phosphoryl oxygen (double-bonded to P) and two ester oxygens (each attached to a carbon that bears a double-bonded oxygen).
Note: This is inherently a heuristic method and it does not capture all stereochemical and connectivity nuances.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

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
    # We look for the C(=O)O substructure in the overall molecule.
    ester_smarts = "C(=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester bond(s); at least 2 are required for a diacyl phosphatidyl moiety"
    
    # ------ Step 2: Identify candidate inositol ring ------
    # Look for a six-membered ring containing only carbons.
    # In a correct 1-phosphatidyl-1D-myo-inositol headgroup, 5 carbons in the ring should have a free -OH,
    # while the 6th carbon carries an oxygen that is substituted by phosphorus.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_inositol = None
    candidate_phospho_O_idx = None  # index of the oxygen substituent on the ring that connects to P
    for ring in ring_info:
        if len(ring) != 6:
            continue  # only interested in six-membered rings

        # Check that all atoms in the ring are carbons.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        
        free_oh_count = 0
        phospho_sub_count = 0
        phospho_O_candidate = None
        # For each ring carbon, look at substituents not part of the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Focus on oxygen substituents.
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen is attached to a phosphorus atom.
                    attached_to_P = any(nbr2.GetAtomicNum() == 15 for nbr2 in nbr.GetNeighbors() if nbr2.GetIdx() != atom.GetIdx())
                    # Check if oxygen carries any hydrogen(s) (i.e. is free hydroxyl).
                    has_H = nbr.GetTotalNumHs() > 0
                    if attached_to_P:
                        phospho_sub_count += 1
                        phospho_O_candidate = nbr
                    elif has_H:
                        free_oh_count += 1
        # For a genuine 1-phosphatidyl-1D-myo-inositol headgroup, we expect exactly 1 phospho substitution and 5 free hydroxyls.
        if phospho_sub_count == 1 and free_oh_count == 5:
            candidate_inositol = ring
            candidate_phospho_O_idx = phospho_O_candidate.GetIdx()
            break

    if candidate_inositol is None:
        return False, "No six-membered inositol ring with 1 phosphate substituent and 5 free hydroxyls found"
    
    # ------ Step 3: Examine connectivity of the phosphate group ------
    # From the inositol oxygen, follow to the phosphorus atom.
    phospho_O = mol.GetAtomWithIdx(candidate_phospho_O_idx)
    phosphorus = None
    for nbr in phospho_O.GetNeighbors():
        if nbr.GetAtomicNum() == 15:
            phosphorus = nbr
            break
    if phosphorus is None:
        return False, "The oxygen substituent on the inositol ring is not attached to any phosphorus atom"
    
    # For a proper phosphatidyl group, phosphorus should have 4 oxygen neighbors.
    p_neighbors = [nbr for nbr in phosphorus.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(p_neighbors) != 4:
        return False, f"Phosphorus atom is bonded to {len(p_neighbors)} oxygen(s); expected 4 (including the inositol oxygen, a phosphoryl O, and 2 ester O's)"

    # Identify the inositol oxygen among the phosphorus neighbors.
    remaining_Os = []
    for o_atom in p_neighbors:
        if o_atom.GetIdx() == candidate_phospho_O_idx:
            continue  # this is the inositol oxygen
        remaining_Os.append(o_atom)
    
    # Now, verify that among the remaining three oxygens one is a phosphoryl oxygen and two are ester oxygens.
    phosphoryl_count = 0
    ester_count = 0
    for o_atom in remaining_Os:
        bond_PO = phosphorus.GetBondBetweenAtoms(phosphorus.GetIdx(), o_atom.GetIdx())
        if bond_PO is None:
            continue
        # Check if the bond is a double bond (P=O).
        if bond_PO.GetBondType() == rdchem.BondType.DOUBLE:
            phosphoryl_count += 1
        else:
            # For a potential ester oxygen, check that this oxygen is connected to a carbon that is itself double-bonded to at least one oxygen.
            is_ester = False
            for nbr in o_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:  # carbon candidate
                    for bond in nbr.GetBonds():
                        # Skip oxygen that is our current o_atom
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                            is_ester = True
                            break
                    if is_ester:
                        break
            if is_ester:
                ester_count += 1
    if phosphoryl_count != 1:
        return False, f"Phosphorus atom is connected to {phosphoryl_count} phosphoryl (P=O) bond(s); expected exactly 1"
    if ester_count != 2:
        return False, f"Phosphorus atom is connected to {ester_count} ester bond(s) (via O-C(=O) groups) instead of 2 expected for a diacyl phosphatidyl moiety"
    
    return True, ("Molecule contains a candidate 1D-myo-inositol headgroup (with 1 phospho substitution on a six-membered carbon ring "
                  "and 5 free hydroxyls) and its phosphate group is connected to a phosphoryl oxygen and 2 ester oxygens as expected.")

# ----- For testing purposes -----
if __name__ == "__main__":
    # Test with one of the provided correct examples.
    test_smiles = ("CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)"
                   "[C@@H](O)[C@H]1O)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC")
    result, explanation = is_1_phosphatidyl_1D_myo_inositol(test_smiles)
    print("Classification:", result)
    print("Explanation:", explanation)