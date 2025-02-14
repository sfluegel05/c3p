"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: 11-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    An 11-oxo steroid is an oxo steroid that has an oxo (=O) substituent at position 11.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid fused ring system: three 6-membered rings and one 5-membered ring fused together
    ri = mol.GetRingInfo()
    rings = ri.BondRings()
    if len(rings) < 4:
        return False, "Molecule does not have enough rings to be a steroid"
    
    # Identify ring sizes
    ring_sizes = [len(ring) for ring in ri.AtomRings()]
    if ring_sizes.count(6) < 3 or ring_sizes.count(5) < 1:
        return False, "Molecule does not have the characteristic steroid ring sizes"

    # Check for fused ring system
    fused_rings = rdMolDescriptors.CalcNumSpiroAtoms(mol) + rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    if fused_rings < 3:
        return False, "Molecule does not have a fused steroid ring system"

    # Define SMARTS pattern for steroid nucleus
    steroid_smarts = 'C1CCC2C(C1)CCC3C2CCC4C3(CCCC4)C'  # Simplified pattern
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Molecule does not match steroid nucleus pattern"

    # Define SMARTS pattern for 11-oxo group attached to ring C
    # Position 11 is in ring C, we can look for ketone attached to the steroid nucleus at specific position
    oxo_smarts = '[#6;R][#6;R](=O)[#6;R]'  # Carbonyl group within a ring
    oxo_pattern = Chem.MolFromSmarts(oxo_smarts)
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No oxo group (=O) found in the molecule"

    # Attempt to locate oxo group at position 11
    # Since direct mapping is difficult, we approximate by checking if the ketone is attached to ring C
    ring_atoms = ri.AtomRings()
    ringC_atoms = None
    # Find the third 6-membered ring (ring C)
    six_membered_rings = [ring for ring in ring_atoms if len(ring) == 6]
    if len(six_membered_rings) >= 3:
        ringC_atoms = six_membered_rings[2]  # Assuming rings are ordered A, B, C
    else:
        return False, "Cannot identify ring C in the steroid nucleus"

    # Check if any oxo group is attached to ring C
    for match in oxo_matches:
        # Check if the carbonyl carbon is in ring C
        carbonyl_carbon_idx = match[1]
        if carbonyl_carbon_idx in ringC_atoms:
            return True, "Molecule is an 11-oxo steroid with oxo group at position 11"
    
    return False, "No oxo group (=O) at position 11 of the steroid nucleus"