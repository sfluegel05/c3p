"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: CHEBI:33913 corrinoid
"""

from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly reduced
    pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking
    alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the corrin nucleus SMARTS pattern
    corrin_smarts = """
    [nH]1ccccc1-cc2c([nH]ccc2-cc3c([nH]ccc3-cc4c[nH]ccc14))    """
    # Remove spaces and newlines
    corrin_smarts = "".join(corrin_smarts.split())
    corrin_pattern = Chem.MolFromSmarts(corrin_smarts)
    if corrin_pattern is None:
        return False, "Invalid SMARTS pattern for corrin nucleus"
    
    # Check if molecule has substructure match with corrin nucleus
    if mol.HasSubstructMatch(corrin_pattern):
        return True, "Molecule contains corrin nucleus"
    else:
        return False, "Molecule does not contain corrin nucleus"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33913',
        'name': 'corrinoid',
        'definition': 'A derivative of the corrin nucleus, which contains four reduced or partly reduced pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.',
        'parents': []
    },
    'config': {
        # Configuration parameters can be included here if necessary
    },
    'message': None,
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Metrics placeholders
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}