"""
Classifies: CHEBI:74222 gamma-lactam
"""
from rdkit import Chem

def is_gamma_lactam(smiles: str):
    """
    Determines if a molecule is a gamma-lactam.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactam, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the gamma-lactam substructure
    gamma_lactam_smarts = '[#6]1-[#6]-[#6]-[#6]-[#7]-1-[#6]=O'
    gamma_lactam = Chem.MolFromSmarts(gamma_lactam_smarts)

    if mol.HasSubstructMatch(gamma_lactam):
        return True, "Contains gamma-lactam substructure"
    else:
        return False, "Does not contain gamma-lactam substructure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:74222',
                          'name': 'gamma-lactam',
                          'definition': 'A lactam in which the amide bond is '
                                        'contained within a five-membered '
                                        'ring, which includes the amide '
                                        'nitrogen and the carbonyl carbon.',
                          'parents': ['CHEBI:24995']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 1-2: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}