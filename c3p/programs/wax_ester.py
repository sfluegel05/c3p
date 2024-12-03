"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ester functional group (COO)
    ester_found = False
    ester_carbon = None
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C' and end_atom.GetTotalDegree() == 3) or \
               (end_atom.GetSymbol() == 'O' and begin_atom.GetSymbol() == 'C' and begin_atom.GetTotalDegree() == 3):
                for neighbor in end_atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetIdx() != begin_atom.GetIdx():
                        ester_found = True
                        ester_carbon = end_atom if end_atom.GetSymbol() == 'C' else begin_atom
                        break
            if ester_found:
                break

    if not ester_found:
        return False, "No ester functional group found"

    # Check for long carbon chains on both sides of the ester group
    def get_chain_length(atom, exclude_idx):
        chain_length = 0
        visited = set()
        stack = [(atom, 0)]
        while stack:
            current_atom, length = stack.pop()
            if current_atom.GetSymbol() == 'C' and current_atom.GetIdx() != exclude_idx:
                chain_length = max(chain_length, length)
                visited.add(current_atom.GetIdx())
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited:
                        stack.append((neighbor, length + 1))
        return chain_length

    if ester_carbon is None:
        return False, "No ester carbon found"

    # Get the two oxygen atoms connected to the ester carbon
    oxygens = [neighbor for neighbor in ester_carbon.GetNeighbors() if neighbor.GetSymbol() == 'O']
    if len(oxygens) != 2:
        return False, "Ester carbon does not have two oxygen neighbors"

    # Check the lengths of the carbon chains attached to the two oxygens
    long_chain_found = False
    for oxygen in oxygens:
        for neighbor in oxygen.GetNeighbors():
            if neighbor.GetIdx() != ester_carbon.GetIdx():
                chain_length = get_chain_length(neighbor, ester_carbon.GetIdx())
                if chain_length >= 8:
                    long_chain_found = True
                    break
        if long_chain_found:
            break

    if not long_chain_found:
        return False, "No long carbon chains found"

    return True, "Molecule is a wax ester"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:10036',
                          'name': 'wax ester',
                          'definition': 'A fatty acid ester resulting from the '
                                        'condensation of the carboxy group of '
                                        'a fatty acid with the alcoholic '
                                        'hydroxy group of a fatty alcohol.',
                          'parents': ['CHEBI:35748', 'CHEBI:73702']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 9,
    'num_false_positives': 1,
    'num_true_negatives': 10,
    'num_false_negatives': 2,
    'precision': 0.9,
    'recall': 0.8181818181818182,
    'f1': 0.8571428571428572,
    'accuracy': None}