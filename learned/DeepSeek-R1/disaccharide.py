"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:15953 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide consists of two monosaccharide units joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Must have exactly two rings (monosaccharide units)
    if num_rings != 2:
        return False, f"Expected 2 rings, found {num_rings}"
    
    # Find bridging oxygen (not in any ring) connecting two different rings
    bridging_oxygens = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            # Check if oxygen is not in any ring
            in_ring = False
            for ring in ring_info.AtomRings():
                if atom.GetIdx() in ring:
                    in_ring = True
                    break
            if in_ring:
                continue
            
            # Check oxygen has exactly two carbon neighbors
            neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(neighbors) != 2:
                continue
            
            # Check neighbors are in different rings
            neighbor_rings = []
            for n in neighbors:
                for ring in ring_info.AtomRings():
                    if n.GetIdx() in ring:
                        neighbor_rings.append(ring)
                        break  # each atom can be in only one main ring
            
            if len(neighbor_rings) == 2 and neighbor_rings[0] != neighbor_rings[1]:
                bridging_oxygens += 1
    
    if bridging_oxygens == 1:
        return True, "Two monosaccharide rings connected by a glycosidic bond"
    else:
        return False, f"Found {bridging_oxygens} bridging oxygens, need exactly 1"