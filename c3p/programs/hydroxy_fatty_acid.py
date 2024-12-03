"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            if any(neighbor.GetAtomicNum() == 8 and neighbor.GetTotalDegree() == 1 for neighbor in atom.GetNeighbors()):
                if any(neighbor.GetAtomicNum() == 8 and neighbor.GetTotalDegree() == 2 for neighbor in atom.GetNeighbors()):
                    carboxylic_acid = True
                    break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for hydroxy groups
    hydroxy_groups = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 1]
    if not hydroxy_groups:
        return False, "No hydroxy groups found"

    # Check for long carbon chain
    carbon_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            carbon_chain_length += 1

    if carbon_chain_length < 4:
        return False, "Carbon chain is too short"

    return True, "Valid hydroxy fatty acid"

# Example usage:
# print(is_hydroxy_fatty_acid("C(C)CCCCCC/C=C\CCCCCCCCCCC(C(=O)O)O"))  # 2-hydroxyerucic acid
# print(is_hydroxy_fatty_acid("C(C(O)=O)C/C=C\CC(/C=C/C=C\C/C=C\C/C=C\C/C=C\CC)O"))  # (4Z,8E,10Z,13Z,16Z,19Z)-7-hydroxydocosahexaenoic acid


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24654',
                          'name': 'hydroxy fatty acid',
                          'definition': 'Any fatty acid carrying one or more '
                                        'hydroxy substituents.',
                          'parents': ['CHEBI:35366', 'CHEBI:35868']},
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
    'num_true_positives': 72,
    'num_false_positives': 17,
    'num_true_negatives': 3,
    'num_false_negatives': 0,
    'precision': 0.8089887640449438,
    'recall': 1.0,
    'f1': 0.8944099378881988,
    'accuracy': None}