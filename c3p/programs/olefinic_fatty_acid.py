"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid (a fatty acid containing at least one C=C double bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group (COOH)
    found_carboxyl = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('C') == 1:
                found_carboxyl = True
                break

    if not found_carboxyl:
        return False, "No carboxylic acid group found"

    # Check for the presence of at least one C=C double bond
    found_double_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'C':
                found_double_bond = True
                break

    if not found_double_bond:
        return False, "No C=C double bond found"

    return True, "Molecule is an olefinic fatty acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:53339',
                          'name': 'olefinic fatty acid',
                          'definition': 'Any fatty acid containing at least '
                                        'one C=C double bond.',
                          'parents': ['CHEBI:27208', 'CHEBI:78840']},
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
    'num_true_positives': 38,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 1,
    'precision': 0.6551724137931034,
    'recall': 0.9743589743589743,
    'f1': 0.7835051546391754,
    'accuracy': None}