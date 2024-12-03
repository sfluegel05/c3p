"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide (Any lipid carrying one or more hydroperoxy substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a hydroperoxy group (-OOH)
    hydroperoxy_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 0:
                    hydroperoxy_found = True
                    break
        if hydroperoxy_found:
            break

    if not hydroperoxy_found:
        return False, "No hydroperoxy group found"

    # Check for the presence of lipid characteristics (long carbon chain)
    carbon_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = atom.GetNeighbors()
            carbon_neighbors = sum(1 for neighbor in neighbors if neighbor.GetSymbol() == 'C')
            if carbon_neighbors >= 2:
                carbon_chain_length += 1

    if carbon_chain_length < 10:  # Arbitrary threshold for a long carbon chain
        return False, "Not a lipid (carbon chain too short)"

    return True, "Lipid hydroperoxide detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61051',
                          'name': 'lipid hydroperoxide',
                          'definition': 'Any lipid carrying one or more '
                                        'hydroperoxy substituents.',
                          'parents': ['CHEBI:18059', 'CHEBI:35924']},
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
    'num_false_positives': 2,
    'num_true_negatives': 12,
    'num_false_negatives': 0,
    'precision': 0.875,
    'recall': 1.0,
    'f1': 0.9333333333333333,
    'accuracy': None}