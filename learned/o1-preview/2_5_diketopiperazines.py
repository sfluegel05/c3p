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

    # SMARTS pattern for 2,5-diketopiperazine core
    pattern = Chem.MolFromSmarts('O=C1NC[C@H](N)C(=O)N1')
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for substructure match
    if mol.HasSubstructMatch(pattern):
        return True, "Contains 2,5-diketopiperazine skeleton"

    # Also check for matches allowing for substitutions and stereochemistry variations
    pattern_generic = Chem.MolFromSmarts('O=C1NCCNC(=O)C1')
    if mol.HasSubstructMatch(pattern_generic):
        return True, "Contains 2,5-diketopiperazine skeleton (generic pattern match)"

    return False, "Does not contain 2,5-diketopiperazine skeleton"


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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Performance metrics can be added here if available
}