"""
Classifies: CHEBI:61689 amino cyclitol
"""
from rdkit import Chem

def is_amino_cyclitol(smiles: str):
    """
    Determines if a molecule is an amino cyclitol (a cyclitol with one or more alcoholic hydroxy groups replaced by substituted or unsubstituted amino groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino cyclitol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify all rings in the molecule
    rings = mol.GetRingInfo()
    if not rings.AtomRings():
        return False, "No rings found in the molecule"

    # Identify all hydroxyl groups and amino groups
    hydroxyl_groups = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors())]
    amino_groups = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']

    if not amino_groups:
        return False, "No amino groups found in the molecule"

    # Check if any hydroxyl group is replaced by an amino group
    for ring in rings.AtomRings():
        ring_atoms = set(ring)
        has_hydroxyl = any(atom_idx in hydroxyl_groups for atom_idx in ring_atoms)
        has_amino = any(atom_idx in amino_groups for atom_idx in ring_atoms)
        if has_hydroxyl and has_amino:
            return True, "Molecule is an amino cyclitol"
    
    return False, "No hydroxyl group replaced by an amino group found in the rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61689',
                          'name': 'amino cyclitol',
                          'definition': 'Any cyclitol having one or more '
                                        'alcoholic hydroxy groups replaced by '
                                        'substituted or unsubstituted amino '
                                        'groups.',
                          'parents': ['CHEBI:23451']},
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
    'num_true_positives': 1,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 11,
    'precision': 1.0,
    'recall': 0.08333333333333333,
    'f1': 0.15384615384615385,
    'accuracy': None}