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

    # Pattern for glycerophosphoserine backbone
    glycerophosphoserine_pattern = Chem.MolFromSmarts("[C@H](CO[P](=O)(O)OC[C@H](N)C(=O)O)O")
    if not mol.HasSubstructMatch(glycerophosphoserine_pattern):
        return False, "No glycerophosphoserine backbone found"

    # Pattern for acyl group at the 1-hydroxy position
    # Ensure the acyl is connected to the oxygen atom at the primary hydroxy position
    acyl_1_hydroxy_pattern = Chem.MolFromSmarts("O[C@@H](COP(O)(=O)O)C(=O)C")
    if not mol.HasSubstructMatch(acyl_1_hydroxy_pattern):
        return False, "Acyl group not correctly found at 1-hydroxy position"

    return True, "Contains 1-acyl-sn-glycero-3-phosphoserine structure"

# Example usage with SMILES strings
smiles_example = "CCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OC[C@H](N)C(O)=O"
result, reason = is_1_acyl_sn_glycero_3_phosphoserine(smiles_example)
print(result, reason)