"""
Classifies: CHEBI:74818 heteroaryl hydroxy compound
"""
from rdkit import Chem

def is_heteroaryl_hydroxy_compound(smiles: str):
    """
    Determines if a molecule is a heteroaryl hydroxy compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a heteroaryl hydroxy compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromaticity and presence of heteroatoms in aromatic rings
    heteroatom_in_ring = False
    hydroxy_group_attached = False

    for ring in mol.GetRingInfo().AtomRings():
        if any(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() not in [6, 1] for atom_idx in ring):
            heteroatom_in_ring = True
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O' and atom.GetIsAromatic():
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'C' and neighbor.GetIsAromatic():
                            hydroxy_group_attached = True
                            break
                if hydroxy_group_attached:
                    break
        if hydroxy_group_attached:
            break

    if not heteroatom_in_ring:
        return False, "No heteroatoms found in aromatic rings"
    
    if not hydroxy_group_attached:
        return False, "No hydroxy groups attached to a heteroarene ring"

    return True, "Molecule is a heteroaryl hydroxy compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:74818',
                          'name': 'heteroaryl hydroxy compound',
                          'definition': 'Any organic aromatic compound having '
                                        'one or more hydroxy groups attached '
                                        'to a heteroarene ring.',
                          'parents': ['CHEBI:33659', 'CHEBI:33822']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 4,
    'num_false_positives': 0,
    'num_true_negatives': 18,
    'num_false_negatives': 14,
    'precision': 1.0,
    'recall': 0.2222222222222222,
    'f1': 0.3636363636363636,
    'accuracy': None}