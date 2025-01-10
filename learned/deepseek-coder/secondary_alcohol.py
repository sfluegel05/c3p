"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: CHEBI:17855 secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom which has two other carbon atoms attached to it.

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

    # Define the substructure pattern for a secondary alcohol
    # The pattern is: [C](-[C])(-[C])-[OH]
    secondary_alcohol_pattern = Chem.MolFromSmarts("[CX4H1][CX4H1][OH]")
    
    # Check if the molecule contains the secondary alcohol pattern
    if mol.HasSubstructMatch(secondary_alcohol_pattern):
        return True, "Contains a hydroxyl group attached to a carbon with two other carbon atoms"
    else:
        return False, "Does not contain a hydroxyl group attached to a carbon with two other carbon atoms"