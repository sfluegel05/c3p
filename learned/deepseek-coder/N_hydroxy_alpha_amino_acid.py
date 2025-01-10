"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
Definition: Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino acid backbone: alpha-carbon with NH2 and COOH
    amino_acid_pattern = Chem.MolFromSmarts("[CX4H]([NH2,NH,NH1,NX3])([CX3](=[OX1])[OX2H1,OX1H0])")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone found"

    # Check for N-hydroxy group (N with at least one OH attached)
    n_hydroxy_pattern = Chem.MolFromSmarts("[NX3][OH]")
    if not mol.HasSubstructMatch(n_hydroxy_pattern):
        return False, "No N-hydroxy group found"

    # Verify the N-hydroxy group is attached to the alpha-carbon
    alpha_carbon_pattern = Chem.MolFromSmarts("[CX4H]([NX3][OH])([CX3](=[OX1])[OX2H1,OX1H0])")
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "N-hydroxy group not attached to alpha-carbon"

    return True, "Contains amino acid backbone with N-hydroxy group attached to alpha-carbon"


__metadata__ = {
    'chemical_class': {
        'id': 'N-hydroxy-alpha-amino-acid',
        'name': 'N-hydroxy-alpha-amino-acid',
        'definition': 'Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.'
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 30,
    'num_false_positives': 2,
    'num_true_negatives': 182407,
    'num_false_negatives': 5,
    'num_negatives': None,
    'precision': 0.9375,
    'recall': 0.8571428571428571,
    'f1': 0.8955223880597015,
    'accuracy': 0.9999615
}