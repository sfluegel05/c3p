"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: CHEBI:21547 N-acetyl-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an N-acyl-amino acid that has acetyl as the acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the acetyl group (CH3-CO-)
    acetyl_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])")
    if not mol.HasSubstructMatch(acetyl_pattern):
        return False, "No acetyl group found"

    # Look for the amino acid backbone (NH2-CHR-COOH or NH-CHR-COOH)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1][CX4][CX3](=[OX1])[OX2H,OX1H0-]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone found"

    # Check if the acetyl group is attached to the nitrogen of the amino acid backbone
    acetyl_nitrogen_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])[NX3;H2,H1][CX4][CX3](=[OX1])[OX2H,OX1H0-]")
    if not mol.HasSubstructMatch(acetyl_nitrogen_pattern):
        return False, "Acetyl group not attached to the nitrogen of the amino acid backbone"

    return True, "Contains an acetyl group attached to the nitrogen of an amino acid backbone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:21547',
                          'name': 'N-acetyl-amino acid',
                          'definition': 'An N-acyl-amino acid that has acetyl as the acyl group.',
                          'parents': ['CHEBI:21547', 'CHEBI:21547']},
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