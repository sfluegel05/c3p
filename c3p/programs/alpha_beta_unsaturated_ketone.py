"""
Classifies: CHEBI:51721 alpha,beta-unsaturated ketone
"""
from rdkit import Chem

def is_alpha_beta_unsaturated_ketone(smiles: str):
    """
    Determines if a molecule is an alpha,beta-unsaturated ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha,beta-unsaturated ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for alpha,beta-unsaturated ketone
    pattern1 = Chem.MolFromSmarts('[#6][C]=[C][C](=O)[#6]')
    pattern2 = Chem.MolFromSmarts('[#6]C#C[C](=O)[#6]')

    if mol.HasSubstructMatch(pattern1):
        return True, "Matches alpha,beta-unsaturated ketone pattern 1"
    elif mol.HasSubstructMatch(pattern2):
        return True, "Matches alpha,beta-unsaturated ketone pattern 2"
    else:
        return False, "Does not match alpha,beta-unsaturated ketone patterns"

# Example usage
smiles = "CC1=CC=CC(=C1)C=CC(=O)C2=CC=C(C=C2)N3CCOCC3"
print(is_alpha_beta_unsaturated_ketone(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51721',
                          'name': 'alpha,beta-unsaturated ketone',
                          'definition': 'A ketone of general formula '
                                        'R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) =/= '
                                        'H) or R(1)C#C-C(=O)R(2) (R(2) =/= H) '
                                        'in which the ketonic C=O function is '
                                        'conjugated to an unsaturated C-C bond '
                                        'at the alpha,beta position.',
                          'parents': ['CHEBI:17087']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Matches alpha,beta-unsaturated ketone pattern 1')\n",
    'num_true_positives': 110,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 9,
    'precision': 0.9821428571428571,
    'recall': 0.9243697478991597,
    'f1': 0.9523809523809523,
    'accuracy': None}