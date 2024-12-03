"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid (PUFA).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PUFA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalValence() == 4:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and 'H' not in neighbors:
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Count the number of double bonds
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bonds += 1

    if double_bonds > 1:
        return True, f"Polyunsaturated fatty acid with {double_bonds} double bonds"
    else:
        return False, "Less than 2 double bonds found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26208',
                          'name': 'polyunsaturated fatty acid',
                          'definition': 'Any fatty acid containing more than '
                                        'one double bond. Acids in this group '
                                        'are reported to have cardioprotective '
                                        'effects; and levels are lowered in '
                                        'chronic fatigue syndrome.',
                          'parents': ['CHEBI:27208']},
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
    'num_true_positives': 75,
    'num_false_positives': 14,
    'num_true_negatives': 6,
    'num_false_negatives': 3,
    'precision': 0.8426966292134831,
    'recall': 0.9615384615384616,
    'f1': 0.8982035928143712,
    'accuracy': None}