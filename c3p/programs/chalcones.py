"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem


def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone (1,3-diphenylpropenone and its derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for chalcones
    chalcone_smarts = 'c1ccccc1C=CC(=O)c2ccccc2'
    pattern = Chem.MolFromSmarts(chalcone_smarts)
    
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule is a chalcone"
    else:
        return False, "Molecule does not match chalcone pattern"



__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23086',
                          'name': 'chalcones',
                          'definition': 'A ketone that is '
                                        '1,3-diphenylpropenone '
                                        '(benzylideneacetophenone), '
                                        'ArCH=CH(=O)Ar, and its derivatives '
                                        'formed by substitution.',
                          'parents': [   'CHEBI:51689',
                                         'CHEBI:72544',
                                         'CHEBI:76224']},
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
    'num_true_positives': 22,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 5,
    'precision': 1.0,
    'recall': 0.8148148148148148,
    'f1': 0.8979591836734693,
    'accuracy': None}