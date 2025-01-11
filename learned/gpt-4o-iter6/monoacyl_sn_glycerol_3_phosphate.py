"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.

    A monoacyl-sn-glycerol 3-phosphate has a glycerol backbone, a phosphate group on carbon 3,
    and a single acyl group at either position 1 or 2.

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

    # Modify the glycerol phosphate pattern to ensure we catch required configurations
    glycerol_phosphate_pattern = Chem.MolFromSmarts("O[C@H](CO)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate correctly configured"

    # Patterns for acyl group at position 1 or 2 of glycerol
    acyl_pos1_pattern = Chem.MolFromSmarts("C(=O)O[C@H](CO)CO")
    acyl_pos2_pattern = Chem.MolFromSmarts("[O-]P(=O)(O)OC[C@H](COP(=O)(O)O)C(=O)")

    pos1_matches = mol.GetSubstructMatches(acyl_pos1_pattern)
    pos2_matches = mol.GetSubstructMatches(acyl_pos2_pattern)

    total_acyl_matches = len(pos1_matches) + len(pos2_matches)
    
    # Ensure there's exactly one acyl group attached
    if total_acyl_matches != 1:
        return False, f"Found {total_acyl_matches} acyl groups, need exactly 1"

    return True, "Contains a glycerol backbone with phosphate and one acyl group attached correctly"

# Example usage:
# result, reason = is_monoacyl_sn_glycerol_3_phosphate("P(OCC(O)COC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)(O)(O)=O")
# print(result, reason)