"""
Classifies: CHEBI:16385 organic sulfide
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (RSR where R is not H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the sulfide group (RSR)
    sulfide_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2:
                # Check if both neighbors are carbon atoms
                if all(neighbor.GetSymbol() == 'C' for neighbor in neighbors):
                    sulfide_found = True
                    break

    if not sulfide_found:
        return False, "No RSR group found"

    # Check if R groups are not hydrogen
    for neighbor in atom.GetNeighbors():
        if any(n.GetSymbol() == 'H' for n in neighbor.GetNeighbors()):
            return False, "One of the R groups is hydrogen"

    return True, "Organic sulfide found"

# Example usage:
# print(is_organic_sulfide("CN1C=CN=C1SC2=C(C=C(C=C2)C=NNC(=O)C3=CC=CC=C3O)[N+](=O)[O-]"))  # Should return True
# print(is_organic_sulfide("CCO"))  # Should return False


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16385',
                          'name': 'organic sulfide',
                          'definition': 'Compounds having the structure RSR (R '
                                        '=/= H). Such compounds were once '
                                        'called thioethers.',
                          'parents': ['CHEBI:26822', 'CHEBI:33261']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 49-50: malformed \\N character escape (<string>, line '
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