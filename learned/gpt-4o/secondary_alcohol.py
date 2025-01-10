"""
Classifies: CHEBI:35681 secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol has a hydroxy group (-OH) attached to a saturated carbon
    atom that is bonded to two other carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a secondary alcohol
    secondary_alcohol_pattern = Chem.MolFromSmarts("[CX4H1]([#6])[#6]O")
    
    # Check for the secondary alcohol pattern
    if mol.HasSubstructMatch(secondary_alcohol_pattern):
        return True, "Contains a hydroxy group attached to a secondary carbon atom"
    
    return False, "Does not meet the criteria for a secondary alcohol"