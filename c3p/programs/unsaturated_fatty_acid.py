"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid (contains at least one C=C or C#C bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0 and atom.GetDegree() == 3:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if 'O' in neighbors and neighbors.count('O') == 2:
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for the presence of at least one C=C or C#C bond
    unsaturated_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            if bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'C':
                unsaturated_bond = True
                break

    if not unsaturated_bond:
        return False, "No C=C or C#C bond found"

    return True, "Contains carboxylic acid group and at least one C=C or C#C bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27208',
                          'name': 'unsaturated fatty acid',
                          'definition': 'Any fatty acid containing at least '
                                        'one C=C or C#C bond.',
                          'parents': ['CHEBI:35366']},
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
    'num_true_positives': 105,
    'num_false_positives': 15,
    'num_true_negatives': 5,
    'num_false_negatives': 2,
    'precision': 0.875,
    'recall': 0.9813084112149533,
    'f1': 0.9251101321585904,
    'accuracy': None}