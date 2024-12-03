"""
Classifies: CHEBI:35681 secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all atoms to find -OH groups attached to a carbon
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 1:  # Hydroxy group
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1 and neighbors[0].GetSymbol() == 'C':  # Attached to a carbon
                carbon = neighbors[0]
                if carbon.GetDegree() == 3:  # Carbon has three bonds (one to OH and two to other carbons)
                    carbon_neighbors = carbon.GetNeighbors()
                    carbon_count = sum(1 for n in carbon_neighbors if n.GetSymbol() == 'C')
                    if carbon_count == 2:
                        return True, "Secondary alcohol found"
    return False, "Not a secondary alcohol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35681',
                          'name': 'secondary alcohol',
                          'definition': 'A secondary alcohol is a compound in '
                                        'which a hydroxy group, -OH, is '
                                        'attached to a saturated carbon atom '
                                        'which has two other carbon atoms '
                                        'attached to it.',
                          'parents': ['CHEBI:30879']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[16:14:25] SMILES Parse Error: syntax error while parsing: '
             'CO[C@@H]1CC(=O)Nc2cc(O)cc(CC\\C=C(C)/[C@H](O)[C@@H](C)[C@H](C\\C=C\\C=C\\C=C\x01)OC(=O)C1(CC1)NC(=O)C1=CCCCC1)c2O\n'
             '[16:14:25] SMILES Parse Error: Failed parsing SMILES '
             "'CO[C@@H]1CC(=O)Nc2cc(O)cc(CC\\C=C(C)/[C@H](O)[C@@H](C)[C@H](C\\C=C\\C=C\\C=C\x01)OC(=O)C1(CC1)NC(=O)C1=CCCCC1)c2O' "
             'for input: '
             "'CO[C@@H]1CC(=O)Nc2cc(O)cc(CC\\C=C(C)/[C@H](O)[C@@H](C)[C@H](C\\C=C\\C=C\\C=C\x01)OC(=O)C1(CC1)NC(=O)C1=CCCCC1)c2O'\n",
    'stdout': '',
    'num_true_positives': 233,
    'num_false_positives': 12,
    'num_true_negatives': 8,
    'num_false_negatives': 1,
    'precision': 0.9510204081632653,
    'recall': 0.9957264957264957,
    'f1': 0.9728601252609604,
    'accuracy': None}