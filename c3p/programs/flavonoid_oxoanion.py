"""
Classifies: CHEBI:60038 flavonoid oxoanion
"""
from rdkit import Chem

def is_flavonoid_oxoanion(smiles: str):
    """
    Determines if a molecule is a flavonoid oxoanion (anion arising from deprotonation of at least one OH group in a flavonoid compound).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid oxoanion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of at least one [O-] group
    has_oxoanion = any(atom.GetSymbol() == 'O' and atom.GetFormalCharge() == -1 for atom in mol.GetAtoms())
    if not has_oxoanion:
        return False, "No oxoanion found"

    # Check for the flavonoid structure (C6-C3-C6 backbone with specific oxygenation pattern)
    # This is a simplified check based on common flavonoid structure patterns
    flavonoid_patterns = [
        'Oc1ccccc1-c2coc3cc(O)ccc3c2=O',  # General flavonoid pattern
        'Oc1ccccc1-c2coc3cc([O-])ccc3c2=O',  # Flavonoid pattern with oxoanion
        'Oc1ccccc1-c2coc3cc([O-])ccc3c2=O',  # Flavonoid pattern with oxoanion
        'Oc1ccccc1-c2coc3cc([O-])ccc3c2=O'   # Flavonoid pattern with oxoanion
    ]

    for pattern in flavonoid_patterns:
        flavonoid_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(flavonoid_pattern):
            return True, "Flavonoid oxoanion"

    return False, "Not a flavonoid structure"

# Example usage:
# smiles = "Oc1cc(ccc1OS([O-])(=O)=O)-c1oc2cc([O-])cc(O)c2c(=O)c1OS([O-])(=O)=O"
# result, reason = is_flavonoid_oxoanion(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60038',
                          'name': 'flavonoid oxoanion',
                          'definition': 'Any anion arising from deprotonation '
                                        'of at least one OH group in a '
                                        'flavonoid compound.',
                          'parents': ['CHEBI:25696']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 1,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 11,
    'precision': 1.0,
    'recall': 0.08333333333333333,
    'f1': 0.15384615384615385,
    'accuracy': None}