"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is composed of two monosaccharide units joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify sugar rings (5 or 6-membered rings with exactly one oxygen)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    sugar_rings = []
    for ring in rings:
        if len(ring) == 5 or len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            o_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            c_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
            if o_count == 1 and c_count == (len(ring) - 1):
                sugar_rings.append(set(ring))
    if len(sugar_rings) != 2:
        return False, f"Found {len(sugar_rings)} monosaccharide rings, need exactly 2"
    
    # Find glycosidic bonds between sugar rings
    glycosidic_bonds = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        a1_idx = atom1.GetIdx()
        a2_idx = atom2.GetIdx()
        # Check for exocyclic oxygen atom connecting two ring carbons from different rings
        if atom1.GetAtomicNum() == 8 and not atom1.IsInRing() and atom2.GetAtomicNum() == 6:
            # Oxygen connected to carbon
            for i, ring1 in enumerate(sugar_rings):
                for j, ring2 in enumerate(sugar_rings):
                    if i >= j:
                        continue
                    if a2_idx in ring1:
                        # Check if oxygen connects to a carbon in the other ring
                        for bond2 in atom1.GetBonds():
                            other_atom = bond2.GetOtherAtom(atom1)
                            if other_atom.GetAtomicNum() == 6 and other_atom.GetIdx() in ring2:
                                glycosidic_bonds.append((atom1.GetIdx(), a2_idx, other_atom.GetIdx()))
        elif atom2.GetAtomicNum() == 8 and not atom2.IsInRing() and atom1.GetAtomicNum() == 6:
            # Oxygen connected to carbon
            for i, ring1 in enumerate(sugar_rings):
                for j, ring2 in enumerate(sugar_rings):
                    if i >= j:
                        continue
                    if a1_idx in ring1:
                        # Check if oxygen connects to a carbon in the other ring
                        for bond2 in atom2.GetBonds():
                            other_atom = bond2.GetOtherAtom(atom2)
                            if other_atom.GetAtomicNum() == 6 and other_atom.GetIdx() in ring2:
                                glycosidic_bonds.append((atom2.GetIdx(), a1_idx, other_atom.GetIdx()))
    if len(glycosidic_bonds) == 0:
        return False, "No glycosidic bond found between monosaccharide units"
    elif len(glycosidic_bonds) > 1:
        return False, f"Found {len(glycosidic_bonds)} glycosidic bonds, need exactly 1"
    
    # Check that the glycosidic bond connects the two sugar rings identified
    glyco_bond = glycosidic_bonds[0]
    atom_indices = set(glyco_bond)
    ring1_indices = sugar_rings[0]
    ring2_indices = sugar_rings[1]
    if (atom_indices & ring1_indices) and (atom_indices & ring2_indices):
        return True, "Contains two monosaccharide units connected by a glycosidic bond"
    else:
        return False, "Glycosidic bond does not connect the two monosaccharide units"

__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'disaccharide',
                              'definition': 'A compound in which two monosaccharides are joined by a glycosidic bond.',
                              'parents': []},
                    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                          'f1_threshold': 0.8,
                          'max_attempts': 5,
                          'max_positive_instances': None,
                          'max_positive_to_test': None,
                          'max_negative_to_test': None,
                          'max_positive_in_prompt': 50,
                          'max_negative_in_prompt': 20,
                          'max_instances_in_prompt': 100,
                          'test_proportion': 0.1},
                    'message': None,
                    'attempt': 1,
                    'success': None,
                    'best': None,
                    'error': '',
                    'stdout': None,
                    'num_true_positives': None,
                    'num_false_positives': None,
                    'num_true_negatives': None,
                    'num_false_negatives': None,
                    'num_negatives': None,
                    'precision': None,
                    'recall': None,
                    'f1': None,
                    'accuracy': None}