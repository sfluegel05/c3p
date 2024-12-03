"""
Classifies: CHEBI:33645 acyclic olefin
"""
from rdkit import Chem

def is_acyclic_olefin(smiles: str):
    """
    Determines if a molecule is an acyclic olefin (acyclic branched or unbranched hydrocarbons having one or more carbon-carbon double bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyclic olefin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is acyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is not acyclic (contains rings)"

    # Check if molecule contains only carbon and hydrogen
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in {'C', 'H'}:
            return False, f"Molecule contains non-hydrocarbon atoms: {atom.GetSymbol()}"

    # Check for the presence of at least one carbon-carbon double bond
    has_double_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            if bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'C':
                has_double_bond = True
                break

    if not has_double_bond:
        return False, "No carbon-carbon double bonds found"

    return True, "Molecule is an acyclic olefin"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33645',
                          'name': 'acyclic olefin',
                          'definition': 'Acyclic branched or unbranched '
                                        'hydrocarbons having one or more '
                                        'carbon-carbon double bond.',
                          'parents': ['CHEBI:33641']},
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
    'num_true_positives': 14,
    'num_false_positives': 0,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}