"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester (any ester of carbamic acid or its N-substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for carbamate ester
    carbamate_ester_smarts = 'O=C(O)N'
    carbamate_ester_pattern = Chem.MolFromSmarts(carbamate_ester_smarts)

    if mol.HasSubstructMatch(carbamate_ester_pattern):
        return True, "Molecule contains carbamate ester group"
    else:
        return False, "Molecule does not contain carbamate ester group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23003',
                          'name': 'carbamate ester',
                          'definition': 'Any ester of carbamic acid or its '
                                        'N-substituted derivatives.',
                          'parents': ['CHEBI:33308']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[20:13:04] SMILES Parse Error: syntax error while parsing: '
             'O=C(CCCCCCCCCC=1NC(/C=C\x02/N=C(C=3NC=CC3)C=C2OC)=CC1)C\n'
             '[20:13:04] SMILES Parse Error: Failed parsing SMILES '
             "'O=C(CCCCCCCCCC=1NC(/C=C\x02/N=C(C=3NC=CC3)C=C2OC)=CC1)C' for "
             'input: '
             "'O=C(CCCCCCCCCC=1NC(/C=C\x02/N=C(C=3NC=CC3)C=C2OC)=CC1)C'\n"
             '[20:13:04] SMILES Parse Error: syntax error while parsing: '
             'O=C(O)/C(/C1=CC=CC=C1)=C\x02/OCO/C2=C(\\C3=CC=CC=C3)/C(=O)NCC(=O)O\n'
             '[20:13:04] SMILES Parse Error: Failed parsing SMILES '
             "'O=C(O)/C(/C1=CC=CC=C1)=C\x02/OCO/C2=C(\\C3=CC=CC=C3)/C(=O)NCC(=O)O' "
             'for input: '
             "'O=C(O)/C(/C1=CC=CC=C1)=C\x02/OCO/C2=C(\\C3=CC=CC=C3)/C(=O)NCC(=O)O'\n",
    'stdout': '',
    'num_true_positives': 29,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}