"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol (C13 to C22).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an alcohol group (-OH)
    alcohol_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
            alcohol_group = True
            break

    if not alcohol_group:
        return False, "No alcohol group found"

    # Count the number of carbon atoms in the longest chain
    carbon_chains = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            chain_length = 1
            visited = set()
            stack = [(atom, chain_length)]
            while stack:
                current_atom, current_length = stack.pop()
                if current_atom.GetIdx() not in visited:
                    visited.add(current_atom.GetIdx())
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in visited:
                            stack.append((neighbor, current_length + 1))
            carbon_chains.append(current_length)

    max_chain_length = max(carbon_chains) if carbon_chains else 0

    if 13 <= max_chain_length <= 22:
        return True, f"Long-chain fatty alcohol with a chain length of {max_chain_length} carbons"
    else:
        return False, f"Chain length of {max_chain_length} carbons is outside the range of C13 to C22"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17135',
                          'name': 'long-chain fatty alcohol',
                          'definition': 'A fatty alcohol with a chain length '
                                        'ranging from C13 to C22.',
                          'parents': ['CHEBI:24026']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[19:48:51] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': '',
    'num_true_positives': 38,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.926829268292683,
    'f1': 0.9620253164556963,
    'accuracy': None}