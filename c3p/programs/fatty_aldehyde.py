"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aldehyde group (carbonyl group with hydrogen)
    aldehyde_group = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2:
                if any(n.GetAtomicNum() == 8 and n.GetTotalDegree() == 1 for n in neighbors):  # Oxygen in carbonyl group
                    if any(n.GetAtomicNum() == 1 for n in neighbors):  # Hydrogen in aldehyde group
                        aldehyde_group = True
                        break
    if not aldehyde_group:
        return False, "No aldehyde group found"

    # Check if the aldehyde group is at one end of the carbon chain
    carbonyl_carbon = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2:
                if any(n.GetAtomicNum() == 8 and n.GetTotalDegree() == 1 for n in neighbors):  # Oxygen in carbonyl group
                    if any(n.GetAtomicNum() == 1 for n in neighbors):  # Hydrogen in aldehyde group
                        carbonyl_carbon = atom
                        break
    if carbonyl_carbon is None:
        return False, "No carbonyl carbon found"

    # Check the rest of the molecule to ensure it's a carbon chain
    chain_atoms = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            chain_atoms.add(atom.GetIdx())
        elif atom.GetAtomicNum() != 1:  # Not hydrogen
            return False, "Non-carbon atoms found in the chain"

    if carbonyl_carbon.GetIdx() not in chain_atoms:
        return False, "Carbonyl carbon not part of the main chain"

    return True, "Valid fatty aldehyde"

# Example usage
smiles_list = [
    "C(=C/C(=O)[H])\\C#CC#CCCCCCC",  # 2-tridecene-4,7-diynal
    "O=CCCCCCCCCCCCCCCCCCCC",  # Eicosanal
    "CCCCCCC=O",  # Octanal
]

for smiles in smiles_list:
    result, reason = is_fatty_aldehyde(smiles)
    print(f"SMILES: {smiles}, Result: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35746',
                          'name': 'fatty aldehyde',
                          'definition': 'An aldehyde formally arising from '
                                        'reduction of the carboxylic acid '
                                        'group of its corresponding fatty '
                                        'acid, having a carbonyl group at one '
                                        'end of the carbon chain.',
                          'parents': ['CHEBI:59768']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: C(=C/C(=O)[H])\\C#CC#CCCCCCC, Result: False, Reason: No '
              'aldehyde group found\n'
              'SMILES: O=CCCCCCCCCCCCCCCCCCCC, Result: False, Reason: No '
              'aldehyde group found\n'
              'SMILES: CCCCCCC=O, Result: False, Reason: No aldehyde group '
              'found\n',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 27,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}