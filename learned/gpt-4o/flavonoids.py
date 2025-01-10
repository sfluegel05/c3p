"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for core flavonoid pattern (C6-C3-C6)
    flavonoid_core_pattern = Chem.MolFromSmarts("c1ccccc1C2=CC=CO2")
    if not mol.HasSubstructMatch(flavonoid_core_pattern):
        return False, "No core C6-C3-C6 pattern found"

    # Check additional oxygen atoms typically present in flavonoids
    oxygens_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygens_count < 2:
        return False, f"Too few oxygen atoms ({oxygens_count}) for a flavonoid"

    # Check for methoxy or hydroxy groups which are common
    methoxy_hydroxy_pattern = Chem.MolFromSmarts("[O;D1]-[CH3]")
    if not mol.HasSubstructMatch(methoxy_hydroxy_pattern):
        return False, "No methoxy or hydroxy groups found"

    return True, "Contains C6-C3-C6 pattern with adequate oxygen substitutions resembling flavonoids"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:CHEBI:47916',
                          'name': 'flavonoid',
                          'definition': 'An organic molecular entity consisting of a polyphenolic compound...',
                          'parents': []},
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