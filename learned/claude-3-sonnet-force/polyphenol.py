"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: CHEBI:36328 polyphenol

Polyphenols are defined as members of the class of phenols that contain 2 or more benzene rings,
each of which is substituted by at least one hydroxy group.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify largest aromatic ring system
    ring_info = mol.GetRingInfo()
    aromatic_rings = [ring_info.AtomRings(i) for i in range(ring_info.NumRings()) if ring_info.IsBondRingAromatic(i)]
    largest_aromatic_ring_system = max(aromatic_rings, key=len, default=[])

    # Check if largest aromatic ring system contains at least 2 benzene rings
    n_benzene_rings = sum(1 for ring in largest_aromatic_ring_system if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring))
    if n_benzene_rings < 2:
        return False, "Less than 2 benzene rings in the largest aromatic ring system"

    # Check for hydroxy groups on benzene rings in the largest aromatic ring system
    hydroxy_on_benzene = False
    for ring in largest_aromatic_ring_system:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                neighbors = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in atom.GetNeighbors()]
                if any(neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1 for neighbor in neighbors):
                    hydroxy_on_benzene = True
                    break
    if not hydroxy_on_benzene:
        return False, "No hydroxy groups on benzene rings in the largest aromatic ring system"

    return True, "Contains 2 or more benzene rings in the largest aromatic ring system, with at least one hydroxy group substituted"