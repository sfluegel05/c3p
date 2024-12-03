"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule is a 6-aminopurine derivative (adenine or substituted adenine).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 6-aminopurine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMILES of 6-aminopurine (adenine)
    adenine_smiles = "Nc1ncnc2[nH]cnc12"
    adenine_mol = Chem.MolFromSmiles(adenine_smiles)

    # Check if the molecule contains the 6-aminopurine substructure
    if mol.HasSubstructMatch(adenine_mol):
        return True, "Molecule contains 6-aminopurine substructure"
    else:
        return False, "Molecule does not contain 6-aminopurine substructure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:20706',
                          'name': '6-aminopurines',
                          'definition': 'Any compound having 6-aminopurine '
                                        '(adenine) as part of its structure.',
                          'parents': ['CHEBI:22527']},
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
    'num_true_positives': 83,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 1,
    'precision': 0.9880952380952381,
    'recall': 0.9880952380952381,
    'f1': 0.9880952380952381,
    'accuracy': None}