"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for two substituent groups attached to the glycerol backbone
    substituent_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() > 1:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O':
                        if any([n.GetSymbol() == 'C' for n in neighbor.GetNeighbors() if n.GetIdx() != atom.GetIdx()]):
                            substituent_count += 1

    if substituent_count < 2:
        return False, "Less than two substituent groups attached to the glycerol backbone"

    return True, "Valid diradylglycerol"

# Example usage
smiles = "O(C(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C[C@@H](O)COC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC"
print(is_diradylglycerol(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76578',
                          'name': 'diradylglycerol',
                          'definition': 'Any lipid that is glycerol bearing '
                                        'two substituent groups - either acyl, '
                                        'alkyl, or alk-1-enyl - at any two of '
                                        'the three possible positions.',
                          'parents': ['CHEBI:35741']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Valid diradylglycerol')\n",
    'num_true_positives': 87,
    'num_false_positives': 9,
    'num_true_negatives': 11,
    'num_false_negatives': 12,
    'precision': 0.90625,
    'recall': 0.8787878787878788,
    'f1': 0.8923076923076922,
    'accuracy': None}