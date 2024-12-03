"""
Classifies: CHEBI:38337 pyrimidone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_pyrimidone(smiles: str):
    """
    Determines if a molecule is a pyrimidone (a pyrimidine carrying one or more oxo substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains a pyrimidine ring
    pyrimidine_pattern = Chem.MolFromSmarts('c1ncncc1')
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine ring found"

    # Check if the molecule contains one or more oxo substituents (C=O)
    oxo_pattern = Chem.MolFromSmarts('C=O')
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No oxo substituents found"

    # Check if the oxo substituents are attached to the pyrimidine ring
    pyrimidine_matches = mol.GetSubstructMatches(pyrimidine_pattern)
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)

    for pyrimidine_match in pyrimidine_matches:
        for oxo_match in oxo_matches:
            if any(atom_idx in pyrimidine_match for atom_idx in oxo_match):
                return True, "Pyrimidone identified"

    return False, "Oxo substituents not attached to pyrimidine ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38337',
                          'name': 'pyrimidone',
                          'definition': 'A pyrimidine carrying one or more oxo '
                                        'substituents.',
                          'parents': ['CHEBI:39447']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 29,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}