"""
Classifies: CHEBI:63567 tetrasaccharide derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrasaccharide_derivative(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains at least four sugar rings
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    sugar_ring_count = 0

    for ring in atom_rings:
        if len(ring) == 5 or len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            if all(atom.GetSymbol() in ['C', 'O'] for atom in ring_atoms):
                oxygen_count = sum(1 for atom in ring_atoms if atom.GetSymbol() == 'O')
                if oxygen_count == 1:
                    sugar_ring_count += 1

    if sugar_ring_count < 4:
        return False, "Less than four sugar rings found"

    # Check for glycosidic bonds (C-O-C)
    glycosidic_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C') or \
               (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'O'):
                glycosidic_bond_count += 1

    if glycosidic_bond_count < 3:
        return False, "Less than three glycosidic bonds found"

    return True, "Tetrasaccharide derivative identified"

# Example usage:
smiles = "[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)NC(C)=O)O[C@H]2[C@@H]([C@H]([C@@H](O[C@@H]2CO)O[C@@H]3[C@H]([C@@H](O[C@@H]([C@@H]3O)CO)O[C@H]4[C@@H]([C@H]([C@@H](O[C@@H]4CO)O)NC(C)=O)O)NC(C)=O)NC(C)=O)O"
result, reason = is_tetrasaccharide_derivative(smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63567',
                          'name': 'tetrasaccharide derivative',
                          'definition': 'An oligosaccharide derivative that is '
                                        'formally obtained from a '
                                        'tetrasaccharide.',
                          'parents': ['CHEBI:63563']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'True Tetrasaccharide derivative identified\n',
    'num_true_positives': 31,
    'num_false_positives': 10,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 0.7560975609756098,
    'recall': 1.0,
    'f1': 0.8611111111111112,
    'accuracy': None}