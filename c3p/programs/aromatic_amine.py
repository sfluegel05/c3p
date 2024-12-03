"""
Classifies: CHEBI:33860 aromatic amine
"""
from rdkit import Chem

def is_aromatic_amine(smiles: str):
    """
    Determines if a molecule is an aromatic amine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    amine_found = False
    aromatic_system_found = False

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetDegree() in [1, 2, 3]:  # Check for primary, secondary, or tertiary amine
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic():
                    amine_found = True
                    aromatic_system_found = True
                    break
            if aromatic_system_found:
                break

    if not amine_found:
        return False, "No amine group found"
    
    if not aromatic_system_found:
        return False, "Amine group not directly linked to an aromatic system"

    return True, "Aromatic amine found"

# Example usage:
# smiles = "Nc1cccc2C(=O)N(Cc12)C1CCC(=O)NC1=O"  # lenalidomide
# print(is_aromatic_amine(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33860',
                          'name': 'aromatic amine',
                          'definition': 'An amino compound in which the amino '
                                        'group is linked directly to an '
                                        'aromatic system.',
                          'parents': ['CHEBI:33659', 'CHEBI:50047']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '[13:49:51] Explicit valence for atom # 10 C, 5, is greater than '
             'permitted\n',
    'stdout': '',
    'num_true_positives': 147,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 25,
    'precision': 0.9865771812080537,
    'recall': 0.8546511627906976,
    'f1': 0.9158878504672896,
    'accuracy': None}