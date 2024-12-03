"""
Classifies: CHEBI:25704 organic sulfate
"""
from rdkit import Chem

def is_organic_sulfate(smiles: str):
    """
    Determines if a molecule is an organic sulfate (compounds of the general formula SO3HOR where R is an organyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    sulfate_pattern = Chem.MolFromSmarts('OS(=O)(=O)O')
    
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate group (SO3H) found"

    # Check if the sulfate group is attached to an organic group
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    for match in sulfate_matches:
        sulfate_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in sulfate_atom.GetNeighbors():
            if neighbor.GetIdx() not in match:
                if neighbor.GetSymbol() == 'C':
                    return True, "Valid organic sulfate"
                else:
                    return False, "Sulfate group not attached to an organic group"

    return False, "Sulfate group not attached to an organic group"

# Example usage
smiles = "OS(=O)(=O)Oc1ccc(cc1)[N+]([O-])=O"
print(is_organic_sulfate(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25704',
                          'name': 'organic sulfate',
                          'definition': 'Compounds of the general formula '
                                        'SO3HOR where R is an organyl group',
                          'parents': ['CHEBI:26820', 'CHEBI:37826']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(False, 'Sulfate group not attached to an organic group')\n",
    'num_true_positives': 83,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 31,
    'precision': 1.0,
    'recall': 0.7280701754385965,
    'f1': 0.8426395939086295,
    'accuracy': None}