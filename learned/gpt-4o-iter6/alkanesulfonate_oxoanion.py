"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is characterized by a sulfonate group directly attached
    to a carbon atom, with the formula pattern R-CS([O-])(=O)=O.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the sulfonate group pattern attached directly to a carbon atom
    sulfonate_pattern = Chem.MolFromSmarts("[CX4,SX4]S(=O)(=O)[O-]")
    
    # Check for the sulfonate group
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "Does not contain the characteristic alkanesulfonate sulfonate group"

    return True, "Contains the characteristic alkanesulfonate sulfonate group"

# Example test cases
# Please include the above example SMILES to test the accuracy and make adjustments as needed.