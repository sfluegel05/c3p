"""
Classifies: CHEBI:76107 deuterated compound
"""
from rdkit import Chem

def is_deuterated_compound(smiles: str):
    """
    Determines if a molecule is a deuterated compound (one or more hydrogen atoms replaced by deuterium).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a deuterated compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of deuterium atoms
    deuterium_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'H' and atom.GetIsotope() == 2]

    if deuterium_atoms:
        return True, f"Deuterated compound with {len(deuterium_atoms)} deuterium atom(s)"
    else:
        return False, "No deuterium atoms found"

# Example usage:
# print(is_deuterated_compound("[2H]C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C(O)=O"))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76107',
                          'name': 'deuterated compound',
                          'definition': 'Any isotopically modified compound '
                                        'that has one or more hydrogen atoms '
                                        'replaced by deuterium.',
                          'parents': ['CHEBI:139358']},
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
    'num_true_positives': 8,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 2,
    'precision': 1.0,
    'recall': 0.8,
    'f1': 0.888888888888889,
    'accuracy': None}