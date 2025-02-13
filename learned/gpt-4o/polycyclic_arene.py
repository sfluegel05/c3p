"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic aromatic hydrocarbon (PAH) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Aromaticity perception
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_PROPERTIES)

    # Get ring information
    ring_info = mol.GetRingInfo()
    
    # Identify aromatic rings
    aromatic_rings = []
    for ring in ring_info.BondRings():
        is_aromatic = True
        for a1, a2 in zip(ring, ring[1:] + ring[:1]):
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if not bond or not bond.GetIsAromatic():
                is_aromatic = False
                break
        if is_aromatic:
            aromatic_rings.append(ring)

    # Check number of aromatic rings
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings, not a polycyclic arene"
    
    # Check for fused aromatic structures (shared edges)
    fused = False
    for i in range(len(aromatic_rings)):
        for j in range(i+1, len(aromatic_rings)):
            # Check if rings share two or more atoms (indicating a fused structure)
            shared_atoms = set(aromatic_rings[i]).intersection(aromatic_rings[j])
            if len(shared_atoms) > 1:
                fused = True
                break
        if fused:
            break

    if not fused:
        return False, "No fused aromatic rings, not a polycyclic arene"
    
    return True, "Structure is consistent with a polycyclic aromatic hydrocarbon"