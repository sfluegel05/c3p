"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycine structure
    glycine_smarts = '[N;H2][CH2][C](=O)[O]'
    glycine = Chem.MolFromSmarts(glycine_smarts)
    if glycine is None:
        return False, "Error in glycine SMARTS pattern"

    # Find glycine substructure in the molecule
    if not mol.HasSubstructMatch(glycine):
        return False, "No glycine substructure found"

    # Define the N-acyl pattern
    n_acyl_smarts = '[C](=O)[N]'
    n_acyl = Chem.MolFromSmarts(n_acyl_smarts)
    if n_acyl is None:
        return False, "Error in N-acyl SMARTS pattern"

    # Find N-acyl substructure in the molecule
    if not mol.HasSubstructMatch(n_acyl):
        return False, "No N-acyl substructure found"

    # Check if N-acyl is attached to glycine
    glycine_matches = mol.GetSubstructMatches(glycine)
    n_acyl_matches = mol.GetSubstructMatches(n_acyl)

    for g_match in glycine_matches:
        for n_match in n_acyl_matches:
            if g_match[0] == n_match[1]:  # N of glycine is the same as N of N-acyl
                return True, "N-acylglycine structure found"

    return False, "N-acyl is not attached to glycine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16180',
                          'name': 'N-acylglycine',
                          'definition': 'An N-acyl-amino acid in which amino '
                                        'acid specified is glycine.',
                          'parents': [   'CHEBI:140325',
                                         'CHEBI:24373',
                                         'CHEBI:51569']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 12,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}