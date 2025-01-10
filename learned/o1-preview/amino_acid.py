"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as a carboxylic acid containing one or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group: C(=O)OH
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O[H]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) == 0:
        return False, "No carboxylic acid group found"
    
    # Look for amino groups: primary or secondary amines (excluding amides)
    amino_group_pattern = Chem.MolFromSmarts("[N;H1,H2;!$(N-C=O)]")
    amino_group_matches = mol.GetSubstructMatches(amino_group_pattern)
    if len(amino_group_matches) == 0:
        return False, "No amino group found"
    
    return True, "Contains carboxylic acid group and one or more amino groups"

__metadata__ = {   'chemical_class': {   'id': None,
                             'name': 'amino acid',
                             'definition': 'A carboxylic acid containing one or more amino groups.',
                             'parents': None},
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
    'attempt': 1,
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