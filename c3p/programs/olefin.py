"""
Classifies: CHEBI:33641 olefin
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_olefin(smiles: str):
    """
    Determines if a molecule is an olefin (acyclic and cyclic hydrocarbons having one or more carbon-carbon double bonds, apart from the formal ones in aromatic compounds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbon-carbon double bonds
    double_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'C':
                double_bonds.append(bond)

    if not double_bonds:
        return False, "No carbon-carbon double bonds found"

    # Check if the double bonds are part of aromatic rings
    aromatic_bonds = []
    for bond in double_bonds:
        if bond.GetIsAromatic():
            aromatic_bonds.append(bond)

    if len(double_bonds) == len(aromatic_bonds):
        return False, "All carbon-carbon double bonds are part of aromatic rings"

    return True, "Contains one or more carbon-carbon double bonds not part of aromatic rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33641',
                          'name': 'olefin',
                          'definition': 'Acyclic and cyclic hydrocarbons '
                                        'having one or more carbon-carbon '
                                        'double bonds, apart from the formal '
                                        'ones in aromatic compounds. The class '
                                        'olefins subsumes alkenes and '
                                        'cycloalkenes and the corresponding '
                                        'polyenes.',
                          'parents': ['CHEBI:24632', 'CHEBI:78840']},
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
    'num_false_positives': 8,
    'num_true_negatives': 12,
    'num_false_negatives': 0,
    'precision': 0.8260869565217391,
    'recall': 1.0,
    'f1': 0.9047619047619047,
    'accuracy': None}