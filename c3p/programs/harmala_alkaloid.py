"""
Classifies: CHEBI:61379 harmala alkaloid
"""
from rdkit import Chem

def is_harmala_alkaloid(smiles: str):
    """
    Determines if a molecule is a harmala alkaloid (based on a 1-methyl-9H-beta-carboline skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a harmala alkaloid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for 1-methyl-9H-beta-carboline skeleton
    harmala_pattern = Chem.MolFromSmarts('C1=CC=C2C(=C1)C3=C(N2)C=CC=N3C')

    if mol.HasSubstructMatch(harmala_pattern):
        return True, "Contains 1-methyl-9H-beta-carboline skeleton"
    else:
        return False, "Does not contain 1-methyl-9H-beta-carboline skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61379',
                          'name': 'harmala alkaloid',
                          'definition': 'Any member of a class of naturally '
                                        'occurring alkaloids based on a '
                                        '1-methyl-9H-beta-carboline skeleton.',
                          'parents': ['CHEBI:22315', 'CHEBI:60834']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[01:49:58] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN(C)[C@@]2(C(=C)[C@]1(CC(C3=C2C4=CC=CC=C4N3)=O)[H])[H]\n'
             '[01:49:58] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN(C)[C@@]2(C(=C)[C@]1(CC(C3=C2C4=CC=CC=C4N3)=O)[H])[H]' "
             'for input: '
             "'C/C=C\x01/CN(C)[C@@]2(C(=C)[C@]1(CC(C3=C2C4=CC=CC=C4N3)=O)[H])[H]'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 81,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}