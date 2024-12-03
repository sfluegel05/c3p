"""
Classifies: CHEBI:37143 organofluorine compound
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound (contains at least one carbon-fluorine bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all bonds in the molecule
    for bond in mol.GetBonds():
        # Check if the bond is between carbon and fluorine
        if (bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'F') or \
           (bond.GetBeginAtom().GetSymbol() == 'F' and bond.GetEndAtom().GetSymbol() == 'C'):
            return True, "Contains at least one carbon-fluorine bond"

    return False, "No carbon-fluorine bond found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37143',
                          'name': 'organofluorine compound',
                          'definition': 'An organofluorine compound is a '
                                        'compound containing at least one '
                                        'carbon-fluorine bond.',
                          'parents': ['CHEBI:17792', 'CHEBI:24062']},
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
    'num_true_positives': 153,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9935064935064936,
    'f1': 0.9967426710097721,
    'accuracy': None}