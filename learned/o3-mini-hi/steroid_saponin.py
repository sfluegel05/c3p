"""
Classifies: CHEBI:61655 steroid saponin
"""
#!/usr/bin/env python
"""
Classifies: Saponin derived from a hydroxysteroid (steroid saponin)
Definition: Any saponin derived from a hydroxysteroid.
This program uses substructure searches to detect a steroid nucleus with at least one hydroxyl group
and the presence of at least one sugar-like ring.
"""

from rdkit import Chem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is defined here as a molecule that contains a hydroxysteroid nucleus (i.e. a steroid
    backbone in which at least one carbon carries an -OH group) covalently linked to a sugar moiety.
    
    The steroid nucleus is approximated using a SMARTS pattern for a fused tetracyclic system (three six-membered rings
    and one five-membered ring; the cyclopentanoperhydrophenanthrene skeleton).
    The sugar is approximated by detecting at least one five- or six-membered ring that contains at least 2 oxygen atoms.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule meets the criteria for a steroid saponin, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a steroid nucleus.
    # This is an approximation of a cyclopentanoperhydrophenanthrene skeleton.
    # Pattern: three fused cyclohexanes and one cyclopentane.
    steroid_core_smarts = "C1CC2CCC3C(C2)CCC4C3(CC1)CC4"
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Failed to construct steroid core pattern"
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid nucleus found"
    
    # Check that the steroid nucleus is hydroxylated.
    # We assume that in at least one match of the steroid core, one of the atoms (typically a carbon)
    # is attached to a hydroxyl (an oxygen having at least one hydrogen neighbor).
    steroid_matches = mol.GetSubstructMatches(steroid_core)
    hydroxyl_on_steroid = False
    for match in steroid_matches:
        # For each atom in the steroid nucleus match, check if any neighbor is an -OH group.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:  # oxygen
                    # Check if this oxygen is bound to a hydrogen
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetAtomicNum() == 1:
                            hydroxyl_on_steroid = True
                            break
                if hydroxyl_on_steroid:
                    break
            if hydroxyl_on_steroid:
                break
        if hydroxyl_on_steroid:
            break
    if not hydroxyl_on_steroid:
        return False, "Steroid nucleus found, but no hydroxyl substituent detected on it"
    
    # Now check for sugar moieties.
    # We define a sugar-like ring as a ring of size 5 or 6 that contains at least 2 oxygen atoms.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    sugar_found = False
    for ring in atom_rings:
        if len(ring) in (5, 6):
            # Count oxygen atoms in the ring
            oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxygen_count >= 2:
                sugar_found = True
                break
    if not sugar_found:
        return False, "No sugar-like ring (glycoside) found"
    
    return True, "Molecule contains a hydroxysteroid nucleus linked to a sugar moiety (steroid saponin)"

# For testing (you can uncomment the following lines)
# test_smiles = "O(C1C(C2(C(C3C(C4(C(=CC3)CC(O)CC4OC5OC(C(OC6OC(C(O)C(O)C6O)CO)C(O)C5O)C)CC2)C1)C)C(C(O)CCC(C)C)C)C7OC(C(O)C(O)C7O)CO"  # Example SMILES from prompt
# result, reason = is_steroid_saponin(test_smiles)
# print(result, reason)