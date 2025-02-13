"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is an amide of a sulfonic acid, characterized by the functional group R-S(=O)₂-NR'₂.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the sulfonamide functional group as a SMARTS pattern
    # Sulfur atom double-bonded to two oxygens and single-bonded to nitrogen: -S(=O)(=O)-N-
    sulfonamide_pattern = Chem.MolFromSmarts('S(=O)(=O)N')
    if sulfonamide_pattern is None:
        return False, "Invalid sulfonamide SMARTS pattern"
    
    # Check if the molecule contains the sulfonamide functional group
    if mol.HasSubstructMatch(sulfonamide_pattern):
        return True, "Contains sulfonamide functional group"
    else:
        return False, "Does not contain sulfonamide functional group"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35358',
        'name': 'sulfonamide',
        'definition': 'An amide of a sulfonic acid RS(=O)₂NR₂.',
        'parents': ['CHEBI:35352']
    }
}