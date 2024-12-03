"""
Classifies: CHEBI:61109 pyrimidine nucleotide-sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_pyrimidine_nucleotide_sugar(smiles: str):
    """
    Determines if a molecule is a pyrimidine nucleotide-sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine nucleotide-sugar, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for pyrimidine nucleobase
    pyrimidine_smarts = Chem.MolFromSmarts('c1cncnc1')
    if not mol.HasSubstructMatch(pyrimidine_smarts):
        return False, "No pyrimidine nucleobase found"

    # Check for nucleotide (phosphate group)
    phosphate_smarts = Chem.MolFromSmarts('P(=O)(O)(O)O')
    if not mol.HasSubstructMatch(phosphate_smarts):
        return False, "No phosphate group found"

    # Check for sugar moiety (furanose or pyranose)
    furanose_smarts = Chem.MolFromSmarts('C1OC(CO)C(O)C1O')
    pyranose_smarts = Chem.MolFromSmarts('C1OC(O)C(O)C(O)C1O')
    if not (mol.HasSubstructMatch(furanose_smarts) or mol.HasSubstructMatch(pyranose_smarts)):
        return False, "No sugar moiety found"

    return True, "Molecule is a pyrimidine nucleotide-sugar"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61109',
                          'name': 'pyrimidine nucleotide-sugar',
                          'definition': 'A nucleotide-sugar whose nucleobase '
                                        'is a pyrimidine.',
                          'parents': ['CHEBI:25609']},
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
    'num_true_positives': 16,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 5,
    'precision': 0.9411764705882353,
    'recall': 0.7619047619047619,
    'f1': 0.8421052631578947,
    'accuracy': None}