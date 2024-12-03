"""
Classifies: CHEBI:38716 carboxylic acid dianion
"""
from rdkit import Chem

def is_carboxylic_acid_dianion(smiles: str):
    """
    Determines if a molecule is a carboxylic acid dianion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxylic acid dianion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate groups (COO-)
    carboxylate_groups = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2:
                oxygens = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']
                if all(o.GetFormalCharge() == -1 for o in oxygens):
                    carboxylate_groups += 1

    # A dianion should have at least two carboxylate groups
    if carboxylate_groups >= 2:
        return True, "Contains at least two carboxylate groups"
    else:
        return False, "Does not contain at least two carboxylate groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38716',
                          'name': 'carboxylic acid dianion',
                          'definition': 'Any dianion containing at least one '
                                        'carboxy group.',
                          'parents': ['CHEBI:29067']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 61-62: malformed \\N character escape (<string>, line '
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