"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern matching for glycerophosphoserine backbone 
    glycerophosphoserine_pattern = Chem.MolFromSmarts("[C@H](CO[P](=O)(O)OC[C@H](N)C(=O)O)O")    
    if not mol.HasSubstructMatch(glycerophosphoserine_pattern):
        return False, "No glycerophosphoserine backbone found"

    # Pattern matching for acyl group at 1-hydroxy position
    acyl_group_pattern = Chem.MolFromSmarts("C(=O)C")
    acyl_matches = mol.GetSubstructMatches(acyl_group_pattern)
    
    # Verify that there is an acyl group connected at the correct position
    glycerol_with_acyl_pattern = Chem.MolFromSmarts("OCC(=O)C")
    if not mol.HasSubstructMatch(glycerol_with_acyl_pattern):
        return False, "Acyl group not found at 1-hydroxy position"
    
    return True, "Contains 1-acyl-sn-glycero-3-phosphoserine structure"

# Example usage, replace with actual SMILES strings
smiles_example = "CCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OC[C@H](N)C(O)=O"
result, reason = is_1_acyl_sn_glycero_3_phosphoserine(smiles_example)
print(result, reason)