"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:47909 2,5-diketopiperazine
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a piperazine-2,5-dione core, which is a six-membered ring
    containing two nitrogen atoms at positions 1 and 4, and two carbonyl groups at positions 2 and 5.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the 2,5-diketopiperazine SMARTS pattern with variable substituents
    pattern = Chem.MolFromSmarts('O=C1NC(=O)CN1')
    
    # Check for substructure match
    if mol.HasSubstructMatch(pattern):
        return True, "Contains 2,5-diketopiperazine core"
    else:
        return False, "Does not contain 2,5-diketopiperazine core"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:47909',
        'name': '2,5-diketopiperazine',
        'definition': 'Any piperazinone that has a piperazine-2,5-dione skeleton.',
        'parents': ['CHEBI:24164', 'CHEBI:48373']
    },
    'config': {
        # Configuration parameters can be added here if needed
    },
    'message': None,
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Performance metrics can be added here if available
}