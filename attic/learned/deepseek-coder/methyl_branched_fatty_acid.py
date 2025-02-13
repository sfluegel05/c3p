"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: CHEBI:76971 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a branched-chain fatty acid containing methyl branches only.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for methyl branches (C attached to exactly 3 other carbons)
    methyl_branch_pattern = Chem.MolFromSmarts("[CX4H3]")
    methyl_branch_matches = mol.GetSubstructMatches(methyl_branch_pattern)
    if len(methyl_branch_matches) == 0:
        return False, "No methyl branches found"

    # Check for long carbon chain (at least 6 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 6:
        return False, "Carbon chain too short to be a fatty acid"

    # Ensure no non-methyl branches (e.g., ethyl, propyl, etc.)
    non_methyl_branch_pattern = Chem.MolFromSmarts("[CX4H2][CX4H2]")
    if mol.HasSubstructMatch(non_methyl_branch_pattern):
        return False, "Non-methyl branches detected"

    return True, "Contains carboxylic acid group, methyl branches, and a long carbon chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76971',
                          'name': 'methyl-branched fatty acid',
                          'definition': 'Any branched-chain fatty acid containing methyl branches only.',
                          'parents': ['CHEBI:76970', 'CHEBI:76972']},
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