"""
Classifies: CHEBI:25248 methyl ester
"""
from rdkit import Chem

def is_methyl_ester(smiles: str):
    """
    Determines if a molecule is a methyl ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the methyl ester group (COC(=O))
    ester_group = Chem.MolFromSmarts("COC(=O)")
    if not mol.HasSubstructMatch(ester_group):
        return False, "No methyl ester group found"

    return True, "Methyl ester group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25248',
                          'name': 'methyl ester',
                          'definition': 'Any carboxylic ester resulting from '
                                        'the formal condensation of a carboxy '
                                        'group with methanol.',
                          'parents': ['CHEBI:33308']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 19-20: malformed \\N character escape (<string>, line '
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