"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all nitrogen atoms in the molecule
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']

    for nitrogen in nitrogen_atoms:
        # Check if nitrogen has exactly one carbon neighbor and two hydrogen atoms
        neighbors = nitrogen.GetNeighbors()
        carbon_neighbors = [neighbor for neighbor in neighbors if neighbor.GetSymbol() == 'C']
        hydrogen_count = sum(1 for neighbor in neighbors if neighbor.GetSymbol() == 'H')

        if len(carbon_neighbors) == 1 and hydrogen_count == 2:
            return True, "Primary amine found"

    return False, "No primary amine found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32877',
                          'name': 'primary amine',
                          'definition': 'A compound formally derived from '
                                        'ammonia by replacing one hydrogen '
                                        'atom by a hydrocarbyl group.',
                          'parents': ['CHEBI:32952', 'CHEBI:50994']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 11-12: malformed \\N character escape (<string>, line '
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