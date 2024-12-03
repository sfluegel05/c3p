"""
Classifies: CHEBI:76224 aromatic ketone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aromatic_ketone(smiles: str):
    """
    Determines if a molecule is an aromatic ketone (a ketone in which the carbonyl group is attached to an aromatic ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of a ketone group
    ketone_group = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDoubleProp('_CIPRank') == 2:  # Oxygen with double bond (carbonyl)
                    ketone_group = True
                    break
            if ketone_group:
                break
    if not ketone_group:
        return False, "No ketone group found"

    # Check if the ketone carbon is attached to an aromatic ring
    aromatic_ketone = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDoubleProp('_CIPRank') == 2:  # Oxygen with double bond (carbonyl)
                    for n in atom.GetNeighbors():
                        if n.GetIsAromatic():
                            aromatic_ketone = True
                            break
            if aromatic_ketone:
                break
    if not aromatic_ketone:
        return False, "Ketone group not attached to an aromatic ring"

    return True, "Aromatic ketone identified"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76224',
                          'name': 'aromatic ketone',
                          'definition': 'A ketone in which the carbonyl group '
                                        'is attached to an aromatic ring.',
                          'parents': ['CHEBI:17087']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': 'key `_CIPRank` exists but does not result in a double value '
             'reason: bad any cast',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}