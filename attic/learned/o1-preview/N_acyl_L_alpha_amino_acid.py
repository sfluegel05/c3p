"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: CHEBI:59949 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid is any L-alpha-amino acid carrying an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for N-acyl-L-alpha-amino acid
    # Pattern represents an alpha carbon connected to:
    # - A carboxyl group C(=O)O
    # - A nitrogen atom which is part of an amide bond (N-C(=O)-R)
    pattern = Chem.MolFromSmarts('[CH1]([N][C](=O)[#6])[C](=O)[O]')

    if mol.HasSubstructMatch(pattern):
        return True, "Molecule is an N-acyl-L-alpha-amino acid"
    else:
        return False, "Molecule is not an N-acyl-L-alpha-amino acid"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:59949',
        'name': 'N-acyl-L-alpha-amino acid',
        'definition': 'Any L-alpha-amino acid carrying an N-acyl substituent.',
        'parents': ['CHEBI:33713', 'CHEBI:83821']
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
    'accuracy': None
}