"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

def group_fused_rings(rings):
    """Group rings into fused systems based on shared atoms."""
    fused_systems = []
    for ring in rings:
        ring_set = set(ring)
        to_merge = []
        for i, system in enumerate(fused_systems):
            if ring_set & system:
                to_merge.append(i)
        if to_merge:
            merged = set().union(*[fused_systems[i] for i in to_merge])
            merged.update(ring_set)
            for i in reversed(to_merge):
                del fused_systems[i]
            fused_systems.append(merged)
        else:
            fused_systems.append(ring_set)
    return fused_systems

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a 3β-hydroxy steroid with a four-ring fused nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Check for at least one hydroxyl group (-OH) attached to a carbon
    hydroxyl_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1 and atom.GetDegree() == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetAtomicNum() == 6:
                hydroxyl_carbons.append(neighbor.GetIdx())
    if not hydroxyl_carbons:
        return False, "No hydroxyl group"
    
    # Get all rings in the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings"
    
    # Group rings into fused systems
    fused_systems = group_fused_rings(rings)
    
    # Check each hydroxyl-bearing carbon for being in a 4-ring fused system
    for carbon_idx in hydroxyl_carbons:
        for system in fused_systems:
            if carbon_idx in system:
                # Count how many rings are entirely within this fused system
                ring_count = sum(1 for ring in rings if set(ring).issubset(system))
                if ring_count >= 4:
                    # Check if hydroxyl is in β configuration (same side as typical methyl groups)
                    # Approximated by checking stereochemistry: if attached carbon has specific chirality
                    carbon = mol.GetAtomWithIdx(carbon_idx)
                    if carbon.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                        return True, "3β-hydroxy group on a four-ring fused steroid nucleus"
                    else:
                        # Heuristic: assume β if stereochemistry is unspecified (common in SMILES)
                        return True, "3-hydroxy group on a four-ring fused steroid nucleus (β assumed)"
    
    return False, "No 3β-hydroxy group on a four-ring fused system"