"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a MUFA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxyl group (-COOH)
    carboxyl_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('C') == 1:
                carboxyl_group = True
                break

    if not carboxyl_group:
        return False, "No carboxyl group found"

    # Check for exactly one double or triple bond in the carbon chain
    double_bonds = 0
    triple_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bonds += 1
        elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            triple_bonds += 1

    if double_bonds + triple_bonds != 1:
        return False, "Molecule does not have exactly one double or triple bond"

    # Check that all other bonds in the carbon chain are single bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() not in [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            return False, "Molecule contains non-single bonds outside of the double/triple bond"

    return True, "Molecule is a monounsaturated fatty acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25413',
                          'name': 'monounsaturated fatty acid',
                          'definition': 'Any fatty acid with one double or '
                                        'triple bond in the fatty acid chain '
                                        'and singly bonded carbon atoms in the '
                                        'rest of the chain. MUFAs have '
                                        'positive effects on the '
                                        'cardiovascular system, and in '
                                        'diabetes treatment.',
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 24,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}