"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    
    # Check if carbon count is between 3 and greater than 27
    if num_carbons < 3:
        return False, f"Number of carbon atoms ({num_carbons}) is less than 3"
    
    # Check if there is an alcohol group (OH)
    has_alcohol = any(atom.GetSymbol() == 'O' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()) for atom in mol.GetAtoms())
    
    if not has_alcohol:
        return False, "No alcohol group found"
    
    # Check if the molecule is aliphatic (no aromatic rings)
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic rings"
    
    return True, "Molecule is a fatty alcohol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24026',
                          'name': 'fatty alcohol',
                          'definition': 'An aliphatic alcohol consisting of a '
                                        'chain of 3 to greater than 27 carbon '
                                        'atoms. Fatty alcohols may be '
                                        'saturated or unsaturated and may be '
                                        'branched or unbranched.',
                          'parents': ['CHEBI:30879', 'CHEBI:61697']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 50,
    'num_false_positives': 17,
    'num_true_negatives': 3,
    'num_false_negatives': 8,
    'precision': 0.746268656716418,
    'recall': 0.8620689655172413,
    'f1': 0.7999999999999999,
    'accuracy': None}