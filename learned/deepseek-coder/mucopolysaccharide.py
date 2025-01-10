"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    Mucopolysaccharides are polysaccharides composed of alternating units of uronic acids and glycosamines,
    and are commonly partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of uronic acid pattern (carboxyl group attached to a sugar)
    uronic_acid_pattern = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[OX2H]")
    if not mol.HasSubstructMatch(uronic_acid_pattern):
        return False, "No uronic acid pattern found"

    # Check for presence of glycosamine pattern (amino group attached to a sugar)
    glycosamine_pattern = Chem.MolFromSmarts("[CX4][NX3H2]")
    if not mol.HasSubstructMatch(glycosamine_pattern):
        return False, "No glycosamine pattern found"

    # Check for sulfuric acid ester pattern (sulfate group attached to a sugar)
    sulfate_pattern = Chem.MolFromSmarts("[OX2][SX4](=[OX1])(=[OX1])[OX2]")
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfuric acid ester pattern found"

    # Check for alternating pattern of uronic acids and glycosamines
    # This is a simplified check and may not catch all cases
    uronic_acid_count = len(mol.GetSubstructMatches(uronic_acid_pattern))
    glycosamine_count = len(mol.GetSubstructMatches(glycosamine_pattern))
    if abs(uronic_acid_count - glycosamine_count) > 1:
        return False, "Uronic acids and glycosamines are not alternating"

    # Check molecular weight - mucopolysaccharides are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for mucopolysaccharide"

    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 10:
        return False, "Too few carbons for mucopolysaccharide"
    if o_count < 5:
        return False, "Too few oxygens for mucopolysaccharide"
    if n_count < 1:
        return False, "Too few nitrogens for mucopolysaccharide"

    return True, "Contains alternating uronic acids and glycosamines with sulfuric acid esterification"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37395',
                          'name': 'mucopolysaccharide',
                          'definition': 'Any of the group of polysaccharides '
                                        'composed of alternating units from '
                                        'uronic acids and glycosamines, and '
                                        'commonly partially esterified with '
                                        'sulfuric acid.',
                          'parents': ['CHEBI:18154', 'CHEBI:18154']},
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