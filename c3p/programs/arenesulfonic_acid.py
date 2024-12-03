"""
Classifies: CHEBI:33555 arenesulfonic acid
"""
from rdkit import Chem

def is_arenesulfonic_acid(smiles: str):
    """
    Determines if a molecule is an arenesulfonic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenesulfonic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of sulfonic acid group (S(=O)(=O)O) and aryl group
    sulfonic_acid_groups = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S' and atom.GetTotalDegree() == 4:
            neighbors = atom.GetNeighbors()
            oxy_count = sum(1 for n in neighbors if n.GetSymbol() == 'O' and n.GetTotalDegree() == 1)
            if oxy_count == 3:
                sulfonic_acid_groups.append(atom)

    if not sulfonic_acid_groups:
        return False, "No sulfonic acid groups found"

    # Check if the sulfonic acid group is attached to an aryl group
    for s_atom in sulfonic_acid_groups:
        for neighbor in s_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetIsAromatic():
                return True, "Arenesulfonic acid found"

    return False, "Sulfonic acid group not attached to an aryl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33555',
                          'name': 'arenesulfonic acid',
                          'definition': 'Organic derivatives of sulfonic acid '
                                        'in which the sulfo group is linked '
                                        'directly to carbon of an aryl group.',
                          'parents': ['CHEBI:33551']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 12-13: malformed \\N character escape (<string>, line '
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