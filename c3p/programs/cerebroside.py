"""
Classifies: CHEBI:23079 cerebroside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_cerebroside(smiles: str):
    """
    Determines if a molecule is a cerebroside (a glycosphingolipid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cerebroside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycosyl part (hexose sugar, e.g., glucose or galactose)
    hexose_smarts = "[C@H]1(O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO)"
    hexose = Chem.MolFromSmarts(hexose_smarts)
    if not mol.HasSubstructMatch(hexose):
        return False, "No hexose sugar found"

    # Check for sphingolipid part (sphingosine backbone)
    sphingosine_smarts = "C[C@@H](O)[C@H](NC=O)CO"
    sphingosine = Chem.MolFromSmarts(sphingosine_smarts)
    if not mol.HasSubstructMatch(sphingosine):
        return False, "No sphingosine backbone found"

    # Check for N-acyl chain
    n_acyl_smarts = "NC(=O)C"
    n_acyl = Chem.MolFromSmarts(n_acyl_smarts)
    if not mol.HasSubstructMatch(n_acyl):
        return False, "No N-acyl chain found"

    return True, "Molecule is classified as a cerebroside"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23079',
                          'name': 'cerebroside',
                          'definition': 'Any member of a group of '
                                        'glycosphingolipids, also known as '
                                        'monoglycosylceramides, which are '
                                        'important components in animal muscle '
                                        'and nerve cell membranes.',
                          'parents': ['CHEBI:17761', 'CHEBI:25513']},
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
    'num_true_positives': 14,
    'num_false_positives': 17,
    'num_true_negatives': 3,
    'num_false_negatives': 10,
    'precision': 0.45161290322580644,
    'recall': 0.5833333333333334,
    'f1': 0.509090909090909,
    'accuracy': None}