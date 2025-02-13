"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule has 6-aminopurine (adenine) as part of its structure.
    A 6-aminopurine is any compound containing an adenine moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains 6-aminopurine, False otherwise
        str: Reason for classification
    """

    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the adenine substructure (6-aminopurine) using SMARTS pattern
    adenine_smarts = 'Nc1ncnc2ncnc12'  # Corrected SMARTS pattern for adenine
    adenine_mol = Chem.MolFromSmarts(adenine_smarts)
    if adenine_mol is None:
        return False, "Error in defining adenine substructure"

    # Check if the molecule contains the adenine substructure
    if mol.HasSubstructMatch(adenine_mol):
        return True, "Contains 6-aminopurine (adenine) substructure"
    else:
        return False, "Does not contain 6-aminopurine (adenine) substructure"


__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': '6-aminopurines',
        'definition': 'Any compound having 6-aminopurine (adenine) as part of its structure.',
        'parents': []
    }
}