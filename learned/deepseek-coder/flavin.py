"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: CHEBI:17621 flavin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    A flavin is a derivative of the dimethylisoalloxazine skeleton with a substituent on the 10 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible SMARTS pattern for the dimethylisoalloxazine core
    flavin_core_pattern = Chem.MolFromSmarts("C1=C2C(=O)N=C3C(=O)NC(=O)C3=NC2=C(C)C=C1")
    if not mol.HasSubstructMatch(flavin_core_pattern):
        return False, "No dimethylisoalloxazine core found"

    # Check for the presence of a substituent at the 10 position (N10)
    # The nitrogen at position 10 should be connected to at least one non-hydrogen atom
    n10_pattern = Chem.MolFromSmarts("[N]1([C]2=C(C)C=C3C(=O)NC(=O)C3=NC2=C1)")
    n10_matches = mol.GetSubstructMatches(n10_pattern)
    if not n10_matches:
        return False, "No nitrogen at position 10 found"

    # Check if the nitrogen at position 10 has a substituent (connected to at least one non-hydrogen atom)
    n10_atom = mol.GetAtomWithIdx(n10_matches[0][0])
    if n10_atom.GetDegree() == 0:
        return False, "Nitrogen at position 10 has no substituent"

    # Check for the presence of at least one substituent on the nitrogen at position 10
    if n10_atom.GetDegree() < 1:
        return False, "Nitrogen at position 10 has no substituent"

    return True, "Contains dimethylisoalloxazine core with a substituent at the 10 position"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17621',
                          'name': 'flavin',
                          'definition': 'A derivative of the '
                                        'dimethylisoalloxazine '
                                        '(7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) '
                                        'skeleton, with a substituent on the 10 position.',
                          'parents': ['CHEBI:24835', 'CHEBI:33281']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}