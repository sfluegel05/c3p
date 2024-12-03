"""
Classifies: CHEBI:26959 thiocarboxylic ester
"""
from rdkit import Chem

def is_thiocarboxylic_ester(smiles: str):
    """
    Determines if a molecule is a thiocarboxylic ester (an ester in which one or both oxygens of an ester group have been replaced by divalent sulfur).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiocarboxylic ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ester_found = False
    thiocarboxylic_ester_found = False

    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'O' and atom2.GetDegree() == 1:
                ester_found = True
                for neighbor in atom1.GetNeighbors():
                    if neighbor.GetSymbol() == 'S' and neighbor.GetDegree() == 2:
                        thiocarboxylic_ester_found = True
                        break
                if thiocarboxylic_ester_found:
                    break
            elif atom1.GetSymbol() == 'O' and atom1.GetDegree() == 1 and atom2.GetSymbol() == 'C':
                ester_found = True
                for neighbor in atom2.GetNeighbors():
                    if neighbor.GetSymbol() == 'S' and neighbor.GetDegree() == 2:
                        thiocarboxylic_ester_found = True
                        break
                if thiocarboxylic_ester_found:
                    break

    if thiocarboxylic_ester_found:
        return True, "Thiocarboxylic ester found"
    elif ester_found:
        return False, "Ester found but no sulfur substitution"
    else:
        return False, "No ester group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26959',
                          'name': 'thiocarboxylic ester',
                          'definition': 'An ester in which one or both oxygens '
                                        'of an ester group have been replaced '
                                        'by divalent sulfur.',
                          'parents': ['CHEBI:33261', 'CHEBI:35701']},
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
    'num_true_positives': 89,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}