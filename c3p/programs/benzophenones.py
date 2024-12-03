"""
Classifies: CHEBI:22726 benzophenones
"""
from rdkit import Chem

def is_benzophenones(smiles: str):
    """
    Determines if a molecule is a benzophenone (aromatic ketone with carbonyl group bonded to 2 phenyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzophenone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all carbonyl groups (C=O)
    carbonyls = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2.0:
                    carbonyls.append(atom)

    if not carbonyls:
        return False, "No carbonyl groups found"

    for carbonyl in carbonyls:
        # Check if the carbonyl carbon is bonded to two aromatic rings
        aromatic_count = 0
        for neighbor in carbonyl.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetIsAromatic():
                # Check if the aromatic carbon is part of a phenyl ring
                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                        if neighbor.GetIdx() in ring:
                            aromatic_count += 1
                            break

        if aromatic_count == 2:
            return True, "Molecule is a benzophenone"

    return False, "Carbonyl group is not bonded to two phenyl groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22726',
                          'name': 'benzophenones',
                          'definition': 'Any aromatic ketone in which the '
                                        'carbonyl group is bonded to 2 phenyl '
                                        'groups.',
                          'parents': ['CHEBI:76224']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 20,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}