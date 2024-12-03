"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin (Members of the class of hydroxyflavan that have a flavan-3-ol skeleton and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for flavan-3-ol skeleton
    flavan_3_ol_skeleton = Chem.MolFromSmarts('Oc1ccc(cc1)[C@H]2Oc3cc(O)ccc3C[C@@H]2O')
    if flavan_3_ol_skeleton is None:
        return None, None

    if not mol.HasSubstructMatch(flavan_3_ol_skeleton):
        return False, "No flavan-3-ol skeleton found"

    return True, "Molecule is a catechin"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23053',
                          'name': 'catechin',
                          'definition': 'Members of the class of hydroxyflavan '
                                        'that have a flavan-3-ol skeleton and '
                                        'its substituted derivatives.',
                          'parents': ['CHEBI:72010']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 7,
    'num_true_negatives': 5,
    'num_false_negatives': 2,
    'precision': 0.5882352941176471,
    'recall': 0.8333333333333334,
    'f1': 0.6896551724137931,
    'accuracy': None}