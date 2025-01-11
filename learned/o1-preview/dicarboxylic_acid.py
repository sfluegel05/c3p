"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:35727 dicarboxylic acid
"""

from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing two carboxy groups (-COOH or -COO-).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Define carboxylic acid group SMARTS pattern
    # Matches both protonated (-COOH) and deprotonated (-COO-) forms
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[O;H1,-1]')
    if carboxylic_acid_pattern is None:
        return False, "Invalid SMARTS pattern for carboxylic acid group"

    # Find matches for carboxylic acid group
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_carboxylic_acid_groups = len(carboxylic_acid_matches)

    if num_carboxylic_acid_groups == 2:
        return True, "Contains exactly two carboxylic acid groups"
    elif num_carboxylic_acid_groups < 2:
        return False, f"Contains only {num_carboxylic_acid_groups} carboxylic acid group(s), less than 2"
    else:
        return False, f"Contains {num_carboxylic_acid_groups} carboxylic acid groups, more than 2"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35727',
                                         'name': 'dicarboxylic acid',
                                         'definition': 'Any carboxylic acid containing two carboxy groups.',
                                         'parents': ['CHEBI:33575']},
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