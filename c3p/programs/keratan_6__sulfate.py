"""
Classifies: CHEBI:18331 keratan 6'-sulfate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_keratan_6__sulfate(smiles: str):
    """
    Determines if a molecule is a keratan 6'-sulfate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a keratan 6'-sulfate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfate group (O=S(=O)(O)O) at the 6'-position
    sulfate_pattern = Chem.MolFromSmarts('OS(=O)(=O)O')
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate group found"

    # Check for keratan structure
    keratan_pattern = Chem.MolFromSmarts('C(CO)O')
    if not mol.HasSubstructMatch(keratan_pattern):
        return False, "No keratan structure found"

    return True, "Molecule is a keratan 6'-sulfate"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18331',
                          'name': "keratan 6'-sulfate",
                          'definition': 'A keratan sulfate with random '
                                        "sulfation at the 6'-position.",
                          'parents': ['CHEBI:60924']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 124-125: malformed \\N character escape (<string>, line '
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