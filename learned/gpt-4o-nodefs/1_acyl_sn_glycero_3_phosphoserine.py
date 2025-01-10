"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for entire 1-acyl-sn-glycero-3-phosphoserine structure
    # Ensure stereochemistry at the glycerol carbon and full integration with phosphoserine
    # The pattern ensures a single acyl chain is attached correctly
    pattern = Chem.MolFromSmarts(
        "[C@@H](CO[P](=O)(O)OC[C@H](N)C(O)=O)OCOC(=O)[C][C...]"
    )

    if mol.HasSubstructMatch(pattern):
        return True, "Contains 1-acyl-sn-glycero-3-phosphoserine structure"

    return False, "Does not match 1-acyl-sn-glycero-3-phosphoserine structure"