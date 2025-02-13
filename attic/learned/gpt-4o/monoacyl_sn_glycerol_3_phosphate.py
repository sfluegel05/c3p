"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    This compound has a glycerol backbone with one acyl group at either position 1 or 2, and a phosphate at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General SMARTS: Glycerol 3-phosphate with a variable position for one acyl group
    # Consider both possible acyl positions: either position 1 or position 2
    sn_glycerol_3_phosphate_pattern = Chem.MolFromSmarts("O[C@H](CO)COP(=O)(O)O")
    acyl_position_1_pattern = Chem.MolFromSmarts("O[C@H](COC(=O))COP(=O)(O)O")
    acyl_position_2_pattern = Chem.MolFromSmarts("OC([C@H](COC(=O)))COP(=O)(O)O")
    
    # Check for the glycerol 3-phosphate backbone
    if not mol.HasSubstructMatch(sn_glycerol_3_phosphate_pattern):
        return False, "No glycerol 3-phosphate backbone found"
    
    # Check for one acyl position attached via ester linkage
    ester_matches_1 = mol.GetSubstructMatches(acyl_position_1_pattern)
    ester_matches_2 = mol.GetSubstructMatches(acyl_position_2_pattern)
    
    # Should have exactly one ester linkage, either in position 1 or position 2, but not both
    total_ester_matches = len(ester_matches_1) + len(ester_matches_2)
    if total_ester_matches != 1:
        return False, f"Found {total_ester_matches} ester linkages, expected exactly 1 attached to glycerol backbone"
    
    return True, "Contains glycerol 3-phosphate backbone and one acyl group esterified to it at the correct position"