"""
Classifies: CHEBI:61655 steroid saponin
"""
#!/usr/bin/env python
"""
Classifies: Saponin derived from a hydroxysteroid (steroid saponin)
Definition: Any saponin derived from a hydroxysteroid.
This improved version uses a more relaxed SMARTS pattern for the steroid nucleus to capture variations
found in natural products, and then ensures that at least one hydroxyl group is directly attached to the nucleus.
It also looks for a sugar-like ring (5- or 6-membered with at least 2 oxygen atoms).
"""

from rdkit import Chem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is defined here as a molecule that contains a hydroxysteroid nucleus (a steroid backbone
    that is hydroxylated) and at least one sugar moiety.
    
    The steroid nucleus is approximated with a relaxed SMARTS pattern for a fused tetracyclic system (3 six-membered rings
    and 1 five-membered ring) typical of many steroids. The hydroxylation is checked by finding at least one -OH group
    attached to an atom of the steroid core.
    
    The sugar moiety is approximated by detecting at least one 5- or 6-membered ring that contains at least 2 oxygen atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria for a steroid saponin, False otherwise.
        str: Reason for the classification.
    """
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more relaxed SMARTS pattern for a steroid nucleus.
    # This pattern aims to cover a cyclopentanoperhydrophenanthrene-like core.
    # It may match molecules with unsaturation or substituents that still represent a steroid skeleton.
    steroid_core_smarts = "C1CCC2C3CCC4C(C3)C2C1C4"
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Failed to construct steroid core pattern"
    
    # Check if the molecule contains a substructure matching the steroid nucleus.
    steroid_matches = mol.GetSubstructMatches(steroid_core)
    if not steroid_matches:
        return False, "No steroid nucleus found"
    
    # Verify that at least one steroid core instance has a hydroxyl (-OH) substituent.
    hydroxyl_on_steroid = False
    for match in steroid_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check neighbors for an oxygen that is also bound to a hydrogen
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:  # oxygen atom
                    # Check if oxygen is bound to at least one hydrogen
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
    
    # Check for at least one sugar-like ring.
    # A sugar-like ring is defined as a 5- or 6-membered ring with at least 2 oxygen atoms.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    sugar_found = False
    for ring in atom_rings:
        if len(ring) in (5, 6):
            oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxygen_count >= 2:
                sugar_found = True
                break
    if not sugar_found:
        return False, "No sugar-like ring (glycoside) found"
    
    return True, "Molecule contains a hydroxysteroid nucleus linked to a sugar moiety (steroid saponin)"

# For testing, you may uncomment lines below to try sample SMILES strings:
# test_smiles = "O(C1C(C2(C(C3C(C4(C(=CC3)CC(O)CC4OC5OC(C(OC6OC(C(O)C(O)C6O)CO)C(O)C5O)C)CC2)C1)C)C(C(O)CCC(C)C)C)C7OC(C(O)C(O)C7O)CO"
# result, reason = is_steroid_saponin(test_smiles)
# print(result, reason)