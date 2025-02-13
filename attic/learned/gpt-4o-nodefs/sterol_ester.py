"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid nucleus pattern (3 hexagonal and 1 pentagonal ring typical of sterols)
    steroid_nucleus_smarts = "C1CC2CCC3C(C1)C2CCC4C3(CCC4)"  # Simplified steroid nucleus pattern
    steroid_nucleus_pattern = Chem.MolFromSmarts(steroid_nucleus_smarts)
    
    if not mol.HasSubstructMatch(steroid_nucleus_pattern):
        return False, "No sterol core structure found"
        
    # Ester group pattern (-C(=O)O-)
    ester_pattern_smarts = "C(=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_pattern_smarts)
    
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Check for long aliphatic chain (multiple CH2)
    long_chain_pattern_smarts = "CCCCCCCCCCCCCCCCCC"  # Pattern for a long aliphatic chain
    long_chain_pattern = Chem.MolFromSmarts(long_chain_pattern_smarts)
    
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long fatty acid chain found"

    return True, "Contains sterol core structure esterified with a long fatty acid chain"

# Example metadata
__metadata__ = {   'chemical_class': {   'id': 'CHEBI:None',
                          'name': 'sterol ester',
                          'definition': 'Chemical entities consisting of a sterol esterified with a fatty acid'},
    'config': {   'llm_model_name': 'rdkit_based_model',
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
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None}