"""
Classifies: CHEBI:22562 anilines
"""
from rdkit import Chem

def is_anilines(smiles: str):
    """
    Determines if a molecule is an aniline (benzene carrying at least one amino substituent and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aniline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of benzene ring
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"

    # Check for the presence of at least one amino group attached to the benzene ring
    amino_pattern = Chem.MolFromSmarts("Nc1ccccc1")
    if mol.HasSubstructMatch(amino_pattern):
        return True, "Contains amino group attached to benzene ring"

    # Check for substituted derivatives
    substituted_amino_pattern = Chem.MolFromSmarts("Nc1ccccc1*")
    if mol.HasSubstructMatch(substituted_amino_pattern):
        return True, "Contains substituted amino group attached to benzene ring"

    return False, "No amino group attached to benzene ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22562',
                          'name': 'anilines',
                          'definition': 'Any  aromatic amine that is benzene '
                                        'carrying at least one amino '
                                        'substituent and its substituted '
                                        'derivatives.',
                          'parents': ['CHEBI:22712', 'CHEBI:33860']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 20-21: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}