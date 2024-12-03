"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester (Any carboxylic ester where the carboxylic acid component is acetic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMILES pattern for acetate ester
    acetate_ester_pattern = Chem.MolFromSmarts("CC(=O)O")

    if mol.HasSubstructMatch(acetate_ester_pattern):
        return True, "Molecule contains the acetate ester group"
    else:
        return False, "Molecule does not contain the acetate ester group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47622',
                          'name': 'acetate ester',
                          'definition': 'Any carboxylic ester where the '
                                        'carboxylic acid component is acetic '
                                        'acid.',
                          'parents': ['CHEBI:33308']},
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
    'num_true_positives': 68,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 0,
    'precision': 0.9315068493150684,
    'recall': 1.0,
    'f1': 0.9645390070921985,
    'accuracy': None}