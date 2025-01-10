"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl compounds
A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    Para-terphenyls have a central benzene ring with two phenyl groups attached in para positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for basic para-terphenyl core structure
    # This pattern matches:
    # - A central 6-membered ring (can be aromatic or not)
    # - Two phenyl groups in para positions
    # - Allows for substitutions
    para_terphenyl_pattern = Chem.MolFromSmarts(
        "[c,C]1[c,C][c,C]([c,C]2[c,C][c,C][c,C][c,C][c,C]2)[c,C][c,C]([c,C]3[c,C][c,C][c,C][c,C][c,C]3)[c,C]1"
    )
    
    if not mol.HasSubstructMatch(para_terphenyl_pattern):
        return False, "No para-terphenyl core structure found"
    
    # Get the match atoms
    match = mol.GetSubstructMatch(para_terphenyl_pattern)
    if not match:
        return False, "Could not map core structure"
    
    # Check that we have three distinct 6-membered rings
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Convert match indices to ring atoms sets for the core structure
    match_atoms = set(match)
    core_rings = []
    for ring in rings:
        ring_set = set(ring)
        if len(ring_set.intersection(match_atoms)) >= 6:  # Ring is part of core
            core_rings.append(ring_set)
    
    # Should find exactly 3 rings that are part of the core
    if len(core_rings) < 3:
        return False, f"Found only {len(core_rings)} rings in core structure"
    
    # Check that the rings are not fused (should share at most 1 atom between any pair)
    for i, ring1 in enumerate(core_rings):
        for j, ring2 in enumerate(core_rings[i+1:], i+1):
            shared = len(ring1.intersection(ring2))
            if shared > 1:
                return False, "Contains fused rings"
    
    # Verify carbons in core rings are sp2 hybridized (aromatic or conjugated)
    for ring in core_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:  # Must be carbon
                continue
            # Check hybridization - must be sp2 or aromatic
            if not (atom.GetIsAromatic() or atom.GetHybridization() == Chem.HybridizationType.SP2):
                return False, "Core structure contains non-sp2 carbons"
    
    # Additional check for connectivity - central ring should connect to both outer rings
    central_ring = None
    for ring in core_rings:
        connections = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring:
                    for other_ring in core_rings:
                        if neighbor.GetIdx() in other_ring:
                            connections += 1
                            break
        if connections >= 2:  # This is likely the central ring
            central_ring = ring
            break
    
    if central_ring is None:
        return False, "Could not identify central ring with correct connectivity"

    return True, "Contains para-terphenyl core structure with appropriate substitution pattern"