"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    
    A monoacyl-sn-glycerol 3-phosphate has a glycerol backbone, a phosphate group on carbon 3,
    and a single acyl group on carbon 1 or 2.

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

    # Refined glycerol backbone with phosphate pattern
    glycerol_phosphate_pattern = Chem.MolFromSmarts("O[C@H]([C@H](COP(=O)(O)O)O)CO")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate correctly configured"

    # Refined acyl group pattern for positions 1 or 2 on the glycerol structure
    acyl_pattern = Chem.MolFromSmarts("C(=O)OC")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    # Verify that exactly one acyl group is correctly attached
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 1"

    return True, "Contains a glycerol backbone with phosphate and one acyl group attached correctly"

# Example usage:
# result, reason = is_monoacyl_sn_glycerol_3_phosphate("P(OCC(O)COC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)(O)(O)=O")
# print(result, reason)