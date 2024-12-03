"""
Classifies: CHEBI:63566 pentasaccharide derivative
"""
from rdkit import Chem

def is_pentasaccharide_derivative(smiles: str):
    """
    Determines if a molecule is a pentasaccharide derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pentasaccharide derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all rings in the molecule
    rings = mol.GetRingInfo().AtomRings()
    
    # Check for the presence of at least 5 saccharide units
    saccharide_count = 0
    for ring in rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if len(ring) == 6 and all(atom.GetSymbol() == 'O' or atom.GetSymbol() == 'C' for atom in atoms):
            saccharide_count += 1

    if saccharide_count >= 5:
        return True, "Contains at least 5 saccharide units"
    
    return False, "Does not contain at least 5 saccharide units"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63566',
                          'name': 'pentasaccharide derivative',
                          'definition': 'An oligosaccharide derivative that is '
                                        'formally obtained from a '
                                        'pentasaccharide.',
                          'parents': ['CHEBI:63563']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 19,
    'num_false_positives': 5,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'precision': 0.7916666666666666,
    'recall': 1.0,
    'f1': 0.8837209302325582,
    'accuracy': None}