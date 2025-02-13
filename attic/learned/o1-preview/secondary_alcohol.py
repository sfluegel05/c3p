"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    
    A secondary alcohol is a compound in which a hydroxy group, -OH, is attached 
    to a saturated carbon atom which has two other carbon atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a secondary alcohol
    # [C;X4;H1;D3] - sp3 carbon with one hydrogen and degree 3
    # Attached to [O;H] - hydroxyl group
    secondary_alcohol_pattern = Chem.MolFromSmarts('[C;X4;H1;D3]([O;H])')

    # Search for matches of the secondary alcohol pattern
    matches = mol.GetSubstructMatches(secondary_alcohol_pattern)
    if matches:
        return True, "Contains a secondary alcohol group"
    else:
        return False, "No secondary alcohol group found"

__metadata__ = {
    'chemical_class': {
        'name': 'secondary alcohol',
        'definition': 'A secondary alcohol is a compound in which a hydroxy group, -OH, '
                      'is attached to a saturated carbon atom which has two other carbon '
                      'atoms attached to it.'
    }
}