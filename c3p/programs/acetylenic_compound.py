"""
Classifies: CHEBI:73474 acetylenic compound
"""
from rdkit import Chem

def is_acetylenic_compound(smiles: str):
    """
    Determines if a molecule is an acetylenic compound (contains a C#C bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetylenic compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all bonds in the molecule
    for bond in mol.GetBonds():
        # Check if the bond is a triple bond between two carbon atoms
        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'C':
                return True, "Contains a C#C bond"

    return False, "No C#C bond found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73474',
                          'name': 'acetylenic compound',
                          'definition': 'Any organic molecule containing a C#C '
                                        'bond.',
                          'parents': ['CHEBI:72695']},
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
    'num_true_positives': 22,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 0,
    'precision': 0.9565217391304348,
    'recall': 1.0,
    'f1': 0.9777777777777777,
    'accuracy': None}