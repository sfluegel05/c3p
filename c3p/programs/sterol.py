"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol (3-hydroxy steroid closely related to cholestan-3-ol).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains a 3-hydroxy group
    has_3_hydroxy = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    idx = neighbor.GetIdx()
                    if any(bond.GetIdx() == 3 for bond in mol.GetBonds() if bond.GetBeginAtomIdx() == idx or bond.GetEndAtomIdx() == idx):
                        has_3_hydroxy = True
                        break
        if has_3_hydroxy:
            break

    if not has_3_hydroxy:
        return False, "No 3-hydroxy group found"

    # Check if the molecule has a steroid backbone (cyclopentanoperhydrophenanthrene)
    steroid_backbone = False
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() >= 4:
        rings = ring_info.AtomRings()
        if any(len(ring) == 5 for ring in rings) and sum(len(ring) == 6 for ring in rings) >= 3:
            steroid_backbone = True

    if not steroid_backbone:
        return False, "No steroid backbone found"

    # Check for additional carbon atoms in the side chain
    additional_carbon_side_chain = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and len(atom.GetNeighbors()) == 4:
            additional_carbon_side_chain = True
            break

    if additional_carbon_side_chain:
        return True, "Sterol with additional carbon atoms in the side chain"
    else:
        return True, "Sterol without additional carbon atoms in the side chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15889',
                          'name': 'sterol',
                          'definition': 'Any 3-hydroxy steroid whose skeleton '
                                        'is closely related to cholestan-3-ol '
                                        '(additional carbon atoms may be '
                                        'present in the side chain).',
                          'parents': ['CHEBI:36834']},
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
    'num_true_positives': 1,
    'num_false_positives': 3,
    'num_true_negatives': 16,
    'num_false_negatives': 18,
    'precision': 0.25,
    'recall': 0.05263157894736842,
    'f1': 0.08695652173913043,
    'accuracy': None}