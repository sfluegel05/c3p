"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is a compound with the general structure RC(=O)H.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldehyde functional group pattern: [CX3H1](=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    
    # Check if the molecule contains the aldehyde pattern
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains the aldehyde functional group (RC(=O)H)"
    else:
        return False, "Does not contain the aldehyde functional group (RC(=O)H)"