"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
"""
Classifies: CHEBI:26607 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid is a carboxylic acid with a straight carbon chain and no double/triple bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for unsaturation (double/triple bonds)
    # Using rdMolDescriptors.CalcNumAliphaticRings to check for unsaturation
    if rdMolDescriptors.CalcNumAliphaticRings(mol) > 0:
        return False, "Molecule contains double or triple bonds (unsaturated)"

    # Check for branching (number of non-terminal carbons with >2 connections)
    branch_pattern = Chem.MolFromSmarts("[CX4;H0,H1,H2]")
    branch_matches = mol.GetSubstructMatches(branch_pattern)
    if len(branch_matches) > 0:
        return False, "Molecule contains branches (not a straight chain)"

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings"

    # Check carbon chain length (at least 4 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short (minimum 4 carbons required)"

    return True, "Straight-chain saturated fatty acid with a carboxylic acid group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26607',
                          'name': 'straight-chain saturated fatty acid',
                          'definition': 'Any saturated fatty acid lacking a side-chain.',
                          'parents': ['CHEBI:26607', 'CHEBI:35366']},
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