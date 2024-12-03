"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax (an organic compound or mixture of compounds that is composed of long-chain molecules and is malleable at ambient temperatures).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for long chains (at least 12 carbon atoms in a row)
    chains = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            chain_length = 0
            visited = set()
            stack = [(atom, 0)]
            while stack:
                current_atom, length = stack.pop()
                if current_atom.GetIdx() in visited:
                    continue
                visited.add(current_atom.GetIdx())
                if current_atom.GetSymbol() == 'C':
                    chain_length = max(chain_length, length + 1)
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'C':
                            stack.append((neighbor, length + 1))
            chains.append(chain_length)
    
    if not any(chain >= 12 for chain in chains):
        return False, "No long carbon chains (at least 12 carbon atoms) found"

    # Check if the molecule has ester functional groups (R-COO-R)
    ester_groups = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'O') or (begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C'):
                if any(neighbor.GetSymbol() == 'C' and neighbor.GetIdx() != begin_atom.GetIdx() and neighbor.GetIdx() != end_atom.GetIdx() for neighbor in begin_atom.GetNeighbors()) and \
                   any(neighbor.GetSymbol() == 'C' and neighbor.GetIdx() != begin_atom.GetIdx() and neighbor.GetIdx() != end_atom.GetIdx() for neighbor in end_atom.GetNeighbors()):
                    ester_groups.append((begin_atom, end_atom))
    
    if not ester_groups:
        return False, "No ester functional groups (R-COO-R) found"

    return True, "Molecule is classified as a wax"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73702',
                          'name': 'wax',
                          'definition': 'A chemical substance that is an '
                                        'organic compound or mixture of '
                                        'compounds that is composed of '
                                        'long-chain molecules and is malleable '
                                        'at ambient temperatures.',
                          'parents': ['CHEBI:59999', 'CHEBI:61697']},
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
    'num_true_positives': 11,
    'num_false_positives': 2,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.8461538461538461,
    'recall': 1.0,
    'f1': 0.9166666666666666,
    'accuracy': None}