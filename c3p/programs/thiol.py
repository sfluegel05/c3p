"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol (contains a -SH group attached to a carbon atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if the atom is sulfur
        if atom.GetSymbol() == 'S':
            # Check if sulfur is bonded to a hydrogen
            neighbors = atom.GetNeighbors()
            if any(neighbor.GetSymbol() == 'H' for neighbor in neighbors):
                # Check if sulfur is bonded to a carbon
                if any(neighbor.GetSymbol() == 'C' for neighbor in neighbors):
                    return True, "Thiol group (-SH) attached to a carbon atom found"
    
    return False, "No thiol group (-SH) attached to a carbon atom found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:29256',
                          'name': 'thiol',
                          'definition': 'An organosulfur compound in which a '
                                        'thiol group, -SH, is attached to a '
                                        'carbon atom of any aliphatic or '
                                        'aromatic moiety.',
                          'parents': ['CHEBI:33261']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 12,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}