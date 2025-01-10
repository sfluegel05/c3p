"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is defined as an unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved enoyl pattern: Capture the double bond specifically between positions 2 and 3
    enoyl_pattern = Chem.MolFromSmarts("[CX3](=O)S[CH2][CH]=[CH].[CX3](=O)S[CH]=[CH2]") 
    
    # More generalized CoA pattern recognizing common features
    coA_pattern = Chem.MolFromSmarts("COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O[C@@H]1O)O)N2C=NC3=C2N=CN=C3N") # Focus on key linkage regions

    # Check for the 2-enoyl double bond pattern
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No 2-enoyl pattern (enoyl thioester linkage with correct positioning) found"

    # Ensure the Coenzyme A (CoA) structure is present
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "No Coenzyme A backbone detected"

    return True, "Contains 2-enoyl feature with Coenzyme A structure"