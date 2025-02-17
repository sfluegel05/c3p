"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: 1-phosphatidyl-1D-myo-inositol
Definition:
    A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer 
    and the phosphatidyl group is located at its position 1.
    
Heuristics in this implementation:
1. Look for at least two ester bonds (simplified as the C(=O)O substructure).
2. Identify a candidate inositol moiety by finding a six-membered ring composed entirely of carbons
   in which the substituents are mostly –OH groups. In genuine 1-phosphatidyl-1D-myo-inositol,
   one of the ring carbons (position 1) carries a phosphatidyl substituent (an oxygen connected to phosphorus)
   and the other 5 carbons have free hydroxyl groups.
   
Note: This is a heuristic approach (using substructure and ring substitution counts) and does not truly
validate stereochemistry. Real stereochemical validation of an inositol isomer is very complex.
"""

from rdkit import Chem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise.
        str: Explanation for the classification result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # --- Criterion 1: Check for diacylglycerol ester bonds ---
    # We expect at least 2 ester bonds (C(=O)O) as part of the phosphatidyl (diacylglycerol) group.
    ester_smarts = "C(=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester bond(s); at least 2 are required for the phosphatidyl moiety"
        
    # --- Criterion 2: Inspect six-membered rings for an inositol headgroup ---
    # In 1-phosphatidyl-1D-myo-inositol the inositol ring should be a cyclohexane with 6 carbons,
    # where 5 of these carbons are substituted with free hydroxyl groups (-OH, i.e. oxygen with at least one hydrogen),
    # and 1 carbon (the 1-position) is substituted with an oxygen that connects to a phosphorus atom.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_found = False
    
    # Loop over all rings in the molecule.
    for ring in ring_info:
        if len(ring) != 6:
            continue  # Only interested in six-membered rings
        # Ensure all atoms in the ring are carbons.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
            
        free_oh_count = 0
        phospho_sub_count = 0
        # For each carbon atom in the ring, inspect its neighbors not in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Get indices of atoms in the ring (for quick checking).
            ring_atom_set = set(ring)
            # Flag to note if this carbon has a substituent that is a phosphate connection.
            has_phospho = False
            # Also count free -OH if present.
            has_free_oh = False
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_atom_set:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen is attached to a phosphorus (phospho substituent)
                    attached_to_P = any(nbr2.GetAtomicNum() == 15 for nbr2 in nbr.GetNeighbors() if nbr2.GetIdx() != atom.GetIdx())
                    # Check if oxygen has (at least one) hydrogen attached.
                    # Note: GetTotalNumHs() returns the total number of (implicit+explicit) hydrogens.
                    has_h = nbr.GetTotalNumHs() > 0
                    if attached_to_P:
                        has_phospho = True
                    elif has_h:
                        has_free_oh = True
            if has_phospho:
                phospho_sub_count += 1
            elif has_free_oh:
                free_oh_count += 1
        # For genuine 1-phosphatidyl-1D-myo-inositol we expect exactly 1 phospho substitution and 5 free -OH.
        if phospho_sub_count == 1 and free_oh_count == 5:
            candidate_found = True
            break  # Found a candidate inositol ring
            
    if not candidate_found:
        return False, "No six-membered inositol ring with 1 phosphate substituent and 5 free hydroxyls found"
        
    # If both criteria are met the molecule is considered to belong to the class.
    return True, "Molecule contains a candidate 1D-myo-inositol headgroup (with 1 phospho and 5 free OH groups) and at least 2 ester bonds consistent with a phosphatidyl moiety"


# Example usage (this part can be removed or commented out when used as a module):
if __name__ == "__main__":
    # A test SMILES string from one of the correct cases
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC"
    result, reason = is_1_phosphatidyl_1D_myo_inositol(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)