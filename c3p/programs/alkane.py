"""
Classifies: CHEBI:18310 alkane
"""
from rdkit import Chem

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane (acyclic branched or unbranched hydrocarbon with formula CnH2n+2).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of only carbon and hydrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in {'C', 'H'}:
            return False, "Contains non-carbon, non-hydrogen atoms"

    # Check for the absence of rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings"

    # Check for the correct number of hydrogens: CnH2n+2
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    num_hydrogens = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    expected_hydrogens = 2 * num_carbons + 2

    if num_hydrogens != expected_hydrogens:
        return False, f"Incorrect number of hydrogens: expected {expected_hydrogens}, found {num_hydrogens}"

    return True, "Valid alkane"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18310',
                          'name': 'alkane',
                          'definition': 'An acyclic branched or unbranched '
                                        'hydrocarbon having the general '
                                        'formula CnH2n+2, and therefore '
                                        'consisting entirely of hydrogen atoms '
                                        'and saturated carbon atoms.',
                          'parents': ['CHEBI:24632', 'CHEBI:33653']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 14,
    'num_false_positives': 0,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}