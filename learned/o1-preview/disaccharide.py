"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: disaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Find rings that could be monosaccharide units (5 or 6 membered rings with oxygen)
    rings = mol.GetRingInfo().AtomRings()
    sugar_rings = []
    for ring in rings:
        if len(ring) == 5 or len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            o_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            if o_count == 1:
                sugar_rings.append(ring)
    if len(sugar_rings) != 2:
        return False, f"Found {len(sugar_rings)} monosaccharide rings, need exactly 2"
    
    # Check for glycosidic bond between the two sugars
    # Glycosidic bond is typically an ether linkage (C-O-C) connecting the rings
    bonds = mol.GetBonds()
    glycosidic_bonds = []
    for bond in bonds:
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and atom1.GetAtomicNum() == 8:
            neighbor1 = [nbr.GetIdx() for nbr in atom1.GetNeighbors() if nbr.GetIdx() != atom2.GetIdx()]
            neighbor2 = [nbr.GetIdx() for nbr in atom2.GetNeighbors() if nbr.GetIdx() != atom1.GetIdx()]
            if atom2.GetAtomicNum() == 6 and atom1.IsInRing() and atom2.IsInRing():
                glycosidic_bonds.append(bond)
        elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE and atom2.GetAtomicNum() == 8:
            neighbor1 = [nbr.GetIdx() for nbr in atom2.GetNeighbors() if nbr.GetIdx() != atom1.GetIdx()]
            neighbor2 = [nbr.GetIdx() for nbr in atom1.GetNeighbors() if nbr.GetIdx() != atom2.GetIdx()]
            if atom1.GetAtomicNum() == 6 and atom2.IsInRing() and atom1.IsInRing():
                glycosidic_bonds.append(bond)
    if len(glycosidic_bonds) == 0:
        return False, "No glycosidic bond found between monosaccharide units"
    
    # Ensure only one glycosidic bond is connecting the two sugars
    if len(glycosidic_bonds) != 1:
        return False, f"Found {len(glycosidic_bonds)} glycosidic bonds, need exactly 1"
    
    # Check that the glycosidic bond connects the two sugar rings identified
    glyco_bond = glycosidic_bonds[0]
    atom1 = glyco_bond.GetBeginAtom()
    atom2 = glyco_bond.GetEndAtom()
    ring1_indices = set(sugar_rings[0])
    ring2_indices = set(sugar_rings[1])
    if (atom1.GetIdx() in ring1_indices and atom2.GetIdx() in ring2_indices) or \
       (atom2.GetIdx() in ring1_indices and atom1.GetIdx() in ring2_indices):
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
        'attempt': 0,
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