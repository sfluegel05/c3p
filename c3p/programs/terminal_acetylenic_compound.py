"""
Classifies: CHEBI:73477 terminal acetylenic compound
"""
from rdkit import Chem

def is_terminal_acetylenic_compound(smiles: str):
    """
    Determines if a molecule is a terminal acetylenic compound (an acetylenic compound where a carbon
    of the C#C moiety is attached to a hydrogen atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a terminal acetylenic compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if the atom is a carbon
        if atom.GetSymbol() == 'C':
            # Check if the carbon is part of a triple bond (C#C)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.TRIPLE:
                    # Check if one of the carbons in the C#C moiety has only one other neighbor (hydrogen)
                    if atom.GetDegree() == 1:
                        return True, "Terminal acetylenic compound found"
                    if neighbor.GetDegree() == 1:
                        return True, "Terminal acetylenic compound found"

    return False, "No terminal acetylenic compound found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73477',
                          'name': 'terminal acetylenic compound',
                          'definition': 'An acetylenic compound which a carbon '
                                        'of the C#C moiety is attached to a '
                                        'hydrogen atom.',
                          'parents': ['CHEBI:73474']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '[02:35:53] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 2,
    'num_true_negatives': 8,
    'num_false_negatives': 0,
    'precision': 0.8333333333333334,
    'recall': 1.0,
    'f1': 0.9090909090909091,
    'accuracy': None}