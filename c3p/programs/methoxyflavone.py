"""
Classifies: CHEBI:25241 methoxyflavone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_methoxyflavone(smiles: str):
    """
    Determines if a molecule is a methoxyflavone (a flavone with at least one methoxy substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methoxyflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the flavone core (2-phenylchromen-4-one)
    flavone_smarts = "O=C1C=CC(=O)C2=C1C=CC=C2"
    flavone_core = Chem.MolFromSmarts(flavone_smarts)
    
    if not mol.HasSubstructMatch(flavone_core):
        return False, "No flavone core found"

    # Check for methoxy groups (OCH3)
    methoxy_smarts = "CO"
    methoxy_group = Chem.MolFromSmarts(methoxy_smarts)
    
    if not mol.HasSubstructMatch(methoxy_group):
        return False, "No methoxy group found"

    return True, "Methoxyflavone found"

# Example usage:
# print(is_methoxyflavone("COc1cc(cc(O)c1O)-c1oc2cc(O)cc(O)c2c(=O)c1O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O"))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25241',
                          'name': 'methoxyflavone',
                          'definition': 'Any member of the class of flavones '
                                        'with at least one methoxy '
                                        'substituent.',
                          'parents': ['CHEBI:24043', 'CHEBI:25698']},
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