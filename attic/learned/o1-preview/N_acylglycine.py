"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: CHEBI:78251 N-acylglycine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is an N-acyl-amino acid where the amino acid is glycine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for N-acylglycine
    # Pattern: Acyl group attached to nitrogen which is connected to CH2-COOH (glycine residue)
    n_acylglycine_pattern = Chem.MolFromSmarts("C(=O)NCC(=O)O")
    if n_acylglycine_pattern is None:
        return False, "Invalid SMARTS pattern for N-acylglycine"

    # Check for substructure match
    if mol.HasSubstructMatch(n_acylglycine_pattern):
        return True, "Contains N-acylglycine substructure"
    else:
        return False, "Does not contain N-acylglycine substructure"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:78251',
        'name': 'N-acylglycine',
        'definition': 'An N-acyl-amino acid in which the amino acid specified is glycine.',
        'parents': ['CHEBI:16549', 'CHEBI:83652']
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
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}