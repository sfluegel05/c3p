"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid (a carboxylic acid containing one or more amino groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    has_carboxylic_acid = False
    has_amino_group = False

    for atom in mol.GetAtoms():
        # Check for carboxylic acid group (COOH or COO-)
        if atom.GetAtomicNum() == 6:  # carbon
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if neighbors.count(8) >= 2:  # two oxygens
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8:  # oxygen
                        if neighbor.GetTotalDegree() == 1:  # OH group
                            has_carboxylic_acid = True
                            break
                if has_carboxylic_acid:
                    break

    for atom in mol.GetAtoms():
        # Check for amino group (NH2, NH, or N)
        if atom.GetAtomicNum() == 7:  # nitrogen
            if any(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors()):  # bonded to hydrogen
                has_amino_group = True
                break

    if has_carboxylic_acid and has_amino_group:
        return True, "Contains both carboxylic acid and amino group"
    elif not has_carboxylic_acid:
        return False, "No carboxylic acid group found"
    elif not has_amino_group:
        return False, "No amino group found"

    return False, "Unknown reason"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33709',
                          'name': 'amino acid',
                          'definition': 'A carboxylic acid containing one or '
                                        'more amino groups.',
                          'parents': ['CHEBI:33575', 'CHEBI:50047']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 10-11: malformed \\N character escape (<string>, line '
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