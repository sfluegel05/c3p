"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid contains an indole skeleton and typically has nitrogen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for indole skeleton pattern (benzene fused to pyrrole)
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]cc2")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole skeleton found"

    # Check for nitrogen atoms (alkaloids typically contain nitrogen)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count == 0:
        return False, "No nitrogen atoms found"

    # Check for basic nitrogen (alkaloids often have basic nitrogen)
    basic_nitrogen = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() > 0:
            basic_nitrogen = True
            break
    if not basic_nitrogen:
        return False, "No basic nitrogen found"

    # Check molecular weight - indole alkaloids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for indole alkaloid"

    return True, "Contains indole skeleton and basic nitrogen atoms"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38958',
                          'name': 'indole alkaloid',
                          'definition': 'An alkaloid containing an indole skeleton.',
                          'parents': ['CHEBI:38958', 'CHEBI:22315']},
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