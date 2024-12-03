"""
Classifies: CHEBI:36054 benzoate ester
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_benzoate_ester(smiles: str):
    """
    Determines if a molecule is a benzoate ester (esters of benzoic acid or substituted benzoic acids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzoate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ester functional group
    ester_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'O':
                for neighbor in atom1.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetIdx() != atom2.GetIdx():
                        ester_found = True
                        break
                if ester_found:
                    break

    if not ester_found:
        return False, "No ester functional group found"

    # Check if the ester is derived from benzoic acid or substituted benzoic acids
    benzoate_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            aromatic_ring = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors() if neighbor.GetIsAromatic()]
            if len(aromatic_ring) >= 6 and all(symbol == 'C' for symbol in aromatic_ring):
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() != atom.GetIdx():
                        for sub_neighbor in neighbor.GetNeighbors():
                            if sub_neighbor.GetSymbol() == 'O' and sub_neighbor.GetIdx() != atom.GetIdx():
                                benzoate_found = True
                                break
                if benzoate_found:
                    break

    if not benzoate_found:
        return False, "No benzoic acid or substituted benzoic acid moiety found"

    return True, "Molecule is a benzoate ester"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36054',
                          'name': 'benzoate ester',
                          'definition': 'Esters of benzoic acid or substituted '
                                        'benzoic acids.',
                          'parents': ['CHEBI:33308', 'CHEBI:62732']},
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
    'num_false_negatives': 83,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}