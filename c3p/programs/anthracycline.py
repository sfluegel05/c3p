"""
Classifies: CHEBI:48120 anthracycline
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthracycline(smiles: str):
    """
    Determines if a molecule is an anthracycline.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthracycline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the tetrahydronaphthacenedione ring structure as SMARTS
    tetrahydronaphthacenedione_smarts = 'O=C1C2=C(O)C3=C(C=CC=C3)C(=O)C4=CC=CC=C4C2=CC1=O'
    tetrahydronaphthacenedione = Chem.MolFromSmarts(tetrahydronaphthacenedione_smarts)
    
    if not mol.HasSubstructMatch(tetrahydronaphthacenedione):
        return False, "No tetrahydronaphthacenedione ring structure found"

    # Define the amino sugar daunosamine as SMARTS
    daunosamine_smarts = 'OC[C@H]1[C@@H](O)[C@H](N)C[C@@H]1O'
    daunosamine = Chem.MolFromSmarts(daunosamine_smarts)
    
    if not mol.HasSubstructMatch(daunosamine):
        return False, "No daunosamine sugar found"

    # Check for glycosidic linkage between tetrahydronaphthacenedione and daunosamine
    matches_tetrahydronaphthacenedione = mol.GetSubstructMatches(tetrahydronaphthacenedione)
    matches_daunosamine = mol.GetSubstructMatches(daunosamine)
    
    for match_tetra in matches_tetrahydronaphthacenedione:
        for match_dauno in matches_daunosamine:
            if any(mol.GetBondBetweenAtoms(match_tetra[i], match_dauno[j]) for i in range(len(match_tetra)) for j in range(len(match_dauno))):
                return True, "Anthracycline structure found with tetrahydronaphthacenedione ring and daunosamine sugar"

    return False, "No glycosidic linkage between tetrahydronaphthacenedione and daunosamine sugar found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48120',
                          'name': 'anthracycline',
                          'definition': 'Anthracyclines are polyketides that '
                                        'have a tetrahydronaphthacenedione '
                                        'ring structure attached by a '
                                        'glycosidic linkage to the amino sugar '
                                        'daunosamine.',
                          'parents': ['CHEBI:26188']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 19,
    'num_false_negatives': 19,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}