"""
Classifies: CHEBI:140345 hydroxy polyunsaturated fatty acid
"""
from rdkit import Chem

def is_hydroxy_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy polyunsaturated fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2:
                carboxylic_acid = True
                break
    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for the presence of hydroxyl groups (OH)
    hydroxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(n.GetSymbol() == 'C' for n in atom.GetNeighbors()))
    if hydroxyl_groups == 0:
        return False, "No hydroxyl groups found"

    # Check for multiple double bonds (polyunsaturation)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bonds < 2:
        return False, "Less than two double bonds found"

    return True, f"Hydroxy polyunsaturated fatty acid with {hydroxyl_groups} hydroxyl groups and {double_bonds} double bonds"

# Example usage
smiles = "C(CCC)C[C@H]([C@@H](C/C=C\\C/C=C\\C/C=C\\CCCC(O)=O)O)O"
print(is_hydroxy_polyunsaturated_fatty_acid(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:140345',
                          'name': 'hydroxy polyunsaturated fatty acid',
                          'definition': 'Any polyunsaturated fatty acid '
                                        'carrying one or more hydroxy '
                                        'substituents.',
                          'parents': ['CHEBI:24654', 'CHEBI:26208']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Hydroxy polyunsaturated fatty acid with 4 hydroxyl "
              "groups and 4 double bonds')\n",
    'num_true_positives': 27,
    'num_false_positives': 11,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.7105263157894737,
    'recall': 1.0,
    'f1': 0.8307692307692308,
    'accuracy': None}