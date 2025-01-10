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

    # Identify sugar rings using SMARTS patterns for pyranose and furanose
    pyranose_smarts = '[C;R]1[C;R][C;R][C;R][C;R][O;R]1'  # 6-membered sugar ring
    furanose_smarts = '[C;R]1[C;R][C;R][C;R][O;R]1'        # 5-membered sugar ring

    pyranose = Chem.MolFromSmarts(pyranose_smarts)
    furanose = Chem.MolFromSmarts(furanose_smarts)

    pyranose_matches = mol.GetSubstructMatches(pyranose)
    furanose_matches = mol.GetSubstructMatches(furanose)

    sugar_ring_atoms = []
    for match in pyranose_matches:
        sugar_ring_atoms.append(set(match))
    for match in furanose_matches:
        sugar_ring_atoms.append(set(match))

    # Ensure there are exactly two sugar rings
    if len(sugar_ring_atoms) != 2:
        return False, f"Found {len(sugar_ring_atoms)} monosaccharide rings, need exactly 2"

    # Identify glycosidic bonds
    # A glycosidic bond is an exocyclic oxygen atom connecting two ring carbons
    glycosidic_bond_smarts = '[C;R][O;!R][C;R]'  # Non-ring oxygen between ring carbons
    glycosidic_bond = Chem.MolFromSmarts(glycosidic_bond_smarts)
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond)

    # Filter glycosidic bonds that connect the two sugar rings
    glycosidic_bonds = []
    for match in glycosidic_bond_matches:
        c1_idx, o_idx, c2_idx = match
        # Check if the two carbons are in different sugar rings
        ring1 = None
        ring2 = None
        for ring in sugar_ring_atoms:
            if c1_idx in ring:
                ring1 = ring
            if c2_idx in ring:
                ring2 = ring
        if ring1 is not None and ring2 is not None and ring1 != ring2:
            glycosidic_bonds.append((c1_idx, o_idx, c2_idx))

    # Ensure there is exactly one glycosidic bond connecting the two sugar rings
    if len(glycosidic_bonds) != 1:
        return False, f"Found {len(glycosidic_bonds)} glycosidic bonds connecting sugar rings, need exactly 1"

    # Check that there are no additional rings (to exclude oligosaccharides)
    total_rings = mol.GetRingInfo().NumRings()
    if total_rings > 2:
        return False, f"Found {total_rings} rings in total, indicates more than two sugar units"

    return True, "Contains two monosaccharide units connected by a glycosidic bond"

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
    'attempt': 2,
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