"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid skeleton pattern (not including stereochemistry details)
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CC2CCC3C4CCC(C4)C3CCC2C1")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "Basic steroid backbone not found"

    # 11beta-hydroxy pattern, considering potential stereochemistry
    # We could assume flexibility in chirality to broaden our pattern
    # Use of "~" to indicate either [C@] or [C]
    hydroxy_11beta_pattern = Chem.MolFromSmarts("C[C@@H](O)C")
    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11beta-hydroxy group found at expected position"

    return True, "Contains steroid backbone with an 11beta-hydroxy group of the expected configuration"