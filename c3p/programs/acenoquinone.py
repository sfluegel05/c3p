"""
Classifies: CHEBI:51285 acenoquinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acenoquinone(smiles: str):
    """
    Determines if a molecule is an acenoquinone (quinones containing an acene fused ring system).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acenoquinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for quinone structure (C=O groups on a ring)
    quinone = False
    for ring in mol.GetRingInfo().AtomRings():
        carbonyl_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2.0:
                        carbonyl_count += 1
        if carbonyl_count >= 2:
            quinone = True
            break

    if not quinone:
        return False, "No quinone structure found"

    # Check for acene (linearly fused aromatic rings)
    acene = False
    ring_info = mol.GetRingInfo()
    aromatic_rings = [ring for ring in ring_info.AtomRings() if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    
    if len(aromatic_rings) >= 3:
        for i in range(len(aromatic_rings) - 2):
            if set(aromatic_rings[i]).intersection(aromatic_rings[i + 1]) and set(aromatic_rings[i + 1]).intersection(aromatic_rings[i + 2]):
                acene = True
                break

    if not acene:
        return False, "No acene fused ring system found"

    return True, "Molecule is an acenoquinone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51285',
                          'name': 'acenoquinone',
                          'definition': 'Quinones containing an acene fused '
                                        'ring system.',
                          'parents': ['CHEBI:36141', 'CHEBI:51269']},
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
    'num_true_positives': 0,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 62,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}