"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone (alpha,beta-unsaturated ketone).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an enone: R1R2C=CR3-C(=O)R4 (R4 != H)
    enone_pattern = Chem.MolFromSmarts('C=CC(=O)[!H]')

    if mol.HasSubstructMatch(enone_pattern):
        return True, "Molecule matches the enone pattern"
    else:
        return False, "Molecule does not match the enone pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51689',
                          'name': 'enone',
                          'definition': 'An alpha,beta-unsaturated ketone of '
                                        'general formula '
                                        'R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) =/= '
                                        'H) in which the C=O function is '
                                        'conjugated to a C=C double bond at '
                                        'the alpha,beta position.',
                          'parents': ['CHEBI:51721', 'CHEBI:78840']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 54-55: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}