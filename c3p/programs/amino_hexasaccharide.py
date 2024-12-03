"""
Classifies: CHEBI:61855 amino hexasaccharide
"""
from rdkit import Chem

def is_amino_hexasaccharide(smiles: str):
    """
    Determines if a molecule is an amino hexasaccharide derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino hexasaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for hexasaccharide backbone
    ring_count = len([ring for ring in mol.GetRingInfo().AtomRings() if len(ring) == 5 or len(ring) == 6])
    if ring_count < 6:
        return False, "The molecule does not have enough rings to be a hexasaccharide"

    # Check for amino groups replacing hydroxy groups
    amino_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
            if 'C' in neighbors:
                amino_found = True
                break

    if not amino_found:
        return False, "No amino group found replacing hydroxy groups"

    return True, "The molecule is an amino hexasaccharide"

# Example usage:
# result, reason = is_amino_hexasaccharide("O[C@@H]1[C@H]([C@@H](O[C@@H]([C@H]1O[C@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O)O[C@]3(O[C@@]([C@@H]([C@H](C3)O)NC(C)=O)([H])[C@@H]([C@@H](CO)O)O)C(O)=O)O)CO[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O[C@H]5[C@@H]([C@H]([C@H]([C@H](O5)CO)O)O[C@H]6[C@@H]([C@H]([C@@H]([C@H](O6)CO)O)O)NC(C)=O)O)O)O)OCCCN)NC(C)=O")
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61855',
                          'name': 'amino hexasaccharide',
                          'definition': 'A hexasaccharide derivative having '
                                        'one or more substituted or '
                                        'unsubstituted amino groups in place '
                                        'of hydroxy groups at unspecified '
                                        'positions.',
                          'parents': ['CHEBI:22483', 'CHEBI:63565']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 3,
    'num_true_negatives': 7,
    'num_false_negatives': 0,
    'precision': 0.7692307692307693,
    'recall': 1.0,
    'f1': 0.8695652173913044,
    'accuracy': None}