"""
Classifies: CHEBI:59268 amino pentasaccharide
"""
from rdkit import Chem

def is_amino_pentasaccharide(smiles: str):
    """
    Determines if a molecule is an amino pentasaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino pentasaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for pentasaccharide structure
    ring_count = 0
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 6:  # Hexose rings
            ring_count += 1
    if ring_count < 5:
        return False, "Not a pentasaccharide (less than 5 hexose rings found)"

    # Check for amino groups
    amino_groups = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
            amino_groups += 1
    if amino_groups == 0:
        return False, "No amino groups found"

    return True, "Valid amino pentasaccharide"

# Example usage:
# smiles = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@@H]3O[C@@H]([C@@H](O)[C@H](OS(O)(=O)=O)[C@H]3O)C(O)=O)[C@H]2O)[C@@H](CO)O[C@H]1O[C@H]1[C@@H](O)[C@@H](CO)O[C@@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@@H]1O"
# print(is_amino_pentasaccharide(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59268',
                          'name': 'amino pentasaccharide',
                          'definition': 'A pentasaccharide derivative having '
                                        'one or more substituted or '
                                        'unsubstituted amino groups in place '
                                        'of hydroxy groups at unspecified '
                                        'positions.',
                          'parents': ['CHEBI:22483', 'CHEBI:63566']},
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
    'num_true_positives': 15,
    'num_false_positives': 0,
    'num_true_negatives': 15,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}