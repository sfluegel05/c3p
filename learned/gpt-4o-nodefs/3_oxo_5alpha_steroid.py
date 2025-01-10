"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for 3-oxo and 5alpha configuration in steroids.
    # 3-Oxo group pattern: C=O at the third carbon position
    # 5alpha configuration with specific cyclopentaphenanthrene core of steroids
    # This pattern captures the carbonyl group at position 3 and required stereochemistry at position 5 
    # considering some flexibility with 5alpha
    oxo_5alpha_pattern = Chem.MolFromSmarts("C1[C@@H]2CC[C@@H]3[C@@H](C(=O))CC[C@@]3(CC[C@@]2([C@@H](C1)C)C1)C1")
    
    if mol.HasSubstructMatch(oxo_5alpha_pattern):
        return True, "Matches the 3-oxo-5alpha-steroid pattern"
    else:
        return False, "Does not contain 3-oxo group in a 5alpha-steroid framework"

    return None, None