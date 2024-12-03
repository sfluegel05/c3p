"""
Classifies: CHEBI:23252 cinnamic acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_cinnamic_acids(smiles: str):
    """
    Determines if a molecule is a cinnamic acid or its substituted derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cinnamic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the cinnamic acid core structure
    cinnamic_acid_core = 'C=CC(=O)O'

    # Check if the molecule contains the cinnamic acid core
    patt = Chem.MolFromSmarts(cinnamic_acid_core)
    if mol.HasSubstructMatch(patt):
        return True, "Contains the cinnamic acid core structure"

    return False, "Does not contain the cinnamic acid core structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23252',
                          'name': 'cinnamic acids',
                          'definition': 'Any alpha,beta-unsaturated '
                                        'monocarboxylic acid based on the '
                                        'cinnamic acid skeleton and its '
                                        'substituted derivatives.',
                          'parents': ['CHEBI:78840', 'CHEBI:79020']},
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
    'num_true_positives': 24,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 9,
    'precision': 0.96,
    'recall': 0.7272727272727273,
    'f1': 0.8275862068965517,
    'accuracy': None}