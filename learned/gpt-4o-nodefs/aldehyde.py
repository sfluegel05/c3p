"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: Aldehyde
"""
from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde contains a carbonyl group (C=O) bonded to a terminal carbon atom.

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
    
    # Define the aldehyde functional group pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")  # This pattern represents a carbon atom with one hydrogen and a double-bonded oxygen.

    # Check if the molecule contains the aldehyde functional group
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains an aldehyde group"
    else:
        return False, "Does not contain an aldehyde group"