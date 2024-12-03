"""
Classifies: CHEBI:76567 polyunsaturated fatty acid anion
"""
from rdkit import Chem

def is_polyunsaturated_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid anion (PUFA anion).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PUFA anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate anion group [C(=O)[O-]]
    carboxylate_anion = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetFormalCharge() == 0:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetFormalCharge() == -1:
                    carboxylate_anion = True
                    break
        if carboxylate_anion:
            break

    if not carboxylate_anion:
        return False, "No carboxylate anion group found"

    # Count the number of C=C double bonds
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'C':
                double_bond_count += 1

    if double_bond_count < 2:
        return False, "Less than two C=C double bonds found"

    return True, f"Polyunsaturated fatty acid anion with {double_bond_count} C=C double bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76567',
                          'name': 'polyunsaturated fatty acid anion',
                          'definition': 'Any unsaturated fatty acid anion '
                                        'containing more than one C-C '
                                        'unsaturated bond.  Major species at '
                                        'pH 7.3.',
                          'parents': ['CHEBI:2580']},
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
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 2,
    'precision': 0.8837209302325582,
    'recall': 0.95,
    'f1': 0.9156626506024096,
    'accuracy': None}