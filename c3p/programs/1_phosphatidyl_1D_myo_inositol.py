"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: 1-phosphatidyl-1D-myo-inositol
Definition: A phosphatidylinositol in which the inositol moiety is the 1D-myo isomer 
            and the phosphatidyl group is located at its position 1.
            
This implementation uses heuristic substructure searches:
1. Looks for a phosphate atom (P with four substituents) that is connected via an O to a six-membered aliphatic ring,
   which is used as a proxy for the inositol headgroup.
2. Searches for at least two ester bonds (C(=O)O) that are expected as part of the diacylglycerol (phosphatidyl) region.
Note: This is an approximative method. Real stereochemical validation is very complex.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Criterion 1: Check for a phosphate group linked to an inositol-like ring ---
    # We require to find at least one phosphorus atom (atomic num 15)
    phosphate_found = False
    inositol_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # phosphorus
            # check that phosphorus has at least one oxygen neighbor that might link to an inositol ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                    # From this oxygen, find a neighbor that is a carbon (and not the phosphorus)
                    for onbr in nbr.GetNeighbors():
                        if onbr.GetIdx() == atom.GetIdx():
                            continue
                        if onbr.GetAtomicNum() == 6:  # carbon atom candidate for inositol headgroup
                            # Retrieve rings for the molecule
                            rings = mol.GetRingInfo().AtomRings()
                            # Check if the candidate carbon is in any ring of size 6
                            for ring in rings:
                                if len(ring) == 6 and onbr.GetIdx() in ring:
                                    # check that all atoms in that ring are carbons 
                                    if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                                        inositol_found = True
                                        break
                            if inositol_found:
                                break
                    if inositol_found:
                        break
            if inositol_found:
                phosphate_found = True
                break  # Found a phosphate with a linkage to a 6-membered carbon ring

    if not phosphate_found or not inositol_found:
        return False, "No phosphate group linked to a six-membered aliphatic (inositol) ring was found"

    # --- Criterion 2: Check for the phosphatidyl (diacylglycerol) moiety ---
    # We expect at least two ester bonds (C(=O)O) present.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Less than 2 ester bonds found (found {len(ester_matches)}), not consistent with a diacylglycerol moiety"
    
    # (Optionally, additional criteria such as checking numbers of rotatable bonds or molecular weight can be added.)
    
    # If both criteria are met, we consider the molecule to belong to the class.
    return True, "Molecule contains a phosphate group linked to a 6-membered carbon (inositol) ring and at least 2 ester bonds consistent with a phosphatidyl moiety"

# Example usage (you can remove or comment these lines when using this as a module):
if __name__ == "__main__":
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC"
    result, reason = is_1_phosphatidyl_1D_myo_inositol(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)