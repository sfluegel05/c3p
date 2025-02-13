"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    This class can have a single acyl group at either position 1 or position 2.

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
    
    # Look for glycerol phosphate backbone without stereo-specific labels
    glycerol_phosphate_pattern = Chem.MolFromSmarts("OC(CO)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate group found"
    
    # Look for presence of a single ester group attached to glycerol backbone
    ester_with_glycerol_pattern = Chem.MolFromSmarts("C(=O)OC(C)COP(=O)(O)O")
    ester_matches_1 = mol.GetSubstructMatches(ester_with_glycerol_pattern)

    # Or at the alternative position
    ester_with_glycerol_pattern_alt = Chem.MolFromSmarts("C(=O)OC(CO)COP(=O)(O)O")
    ester_matches_2 = mol.GetSubstructMatches(ester_with_glycerol_pattern_alt)
    
    # Both patterns should not match at once, and at least one match should be found
    if len(ester_matches_1) > 0 and len(ester_matches_2) > 0:
        return False, "Multiple ester attachments found"
    elif len(ester_matches_1) == 0 and len(ester_matches_2) == 0:
        return False, "No ester linkage typical of an acyl group found"
    
    # Ensure only one acyl group is present
    acyl_group_count = len(ester_matches_1) + len(ester_matches_2)
    if acyl_group_count != 1:
        return False, f"Found {acyl_group_count} acyl groups, need exactly 1"

    return True, "Contains glycerol backbone with one acyl group and a phosphate group"