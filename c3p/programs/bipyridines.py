"""
Classifies: CHEBI:50511 bipyridines
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_bipyridines(smiles: str):
    """
    Determines if a molecule is a bipyridine (compounds containing a bipyridine group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bipyridine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define bipyridine substructure
    bipyridine_smarts = 'n1ccccc1-n2ccccc2'
    bipyridine_mol = Chem.MolFromSmarts(bipyridine_smarts)

    if bipyridine_mol is None:
        return False, "Error in defining bipyridine substructure"

    # Check if the molecule contains the bipyridine substructure
    if mol.HasSubstructMatch(bipyridine_mol):
        return True, "Contains bipyridine group"
    
    # Check if the molecule contains two separate pyridine rings
    pyridine_smarts = 'n1ccccc1'
    pyridine_mol = Chem.MolFromSmarts(pyridine_smarts)
    pyridine_matches = mol.GetSubstructMatches(pyridine_mol)

    if len(pyridine_matches) >= 2:
        return True, "Contains two pyridine rings, possible bipyridine group"

    return False, "Does not contain bipyridine group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50511',
                          'name': 'bipyridines',
                          'definition': 'Compounds containing a bipyridine '
                                        'group.',
                          'parents': [   'CHEBI:36820',
                                         'CHEBI:38101',
                                         'CHEBI:64459']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 14,
    'num_false_positives': 2,
    'num_true_negatives': 13,
    'num_false_negatives': 1,
    'precision': 0.875,
    'recall': 0.9333333333333333,
    'f1': 0.9032258064516129,
    'accuracy': None}