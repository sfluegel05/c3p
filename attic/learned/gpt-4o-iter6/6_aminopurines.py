"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule belongs to the 6-aminopurines class based on its SMILES string.
    6-aminopurine (adenine) is characterized by a specific purine structure with an amino group at the 6th position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains 6-aminopurine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved adenine (6-aminopurine) substructure pattern
    adenine_pattern = Chem.MolFromSmarts("c1ncnc2c1ncn2")
    
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Does not contain 6-aminopurine (adenine) structure"
    
    return True, "Contains 6-aminopurine (adenine) as part of its structure"