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
    try:
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    except Exception as e:
        return False, f"Failed to sanitize and process aromatics: {str(e)}"

    # Get ring information
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized():
        return False, "Ring information is not initialized"
    
    # Identify aromatic rings
    aromatic_rings = [
        ring for ring in ring_info.AtomRings() 
        if all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring)
    ]

    # Check number of rings
    if len(aromatic_rings) < 2:
        return False, f"Less than two aromatic rings, not a polycyclic arene"
    
    # Check for fused aromatic structures (shared bonds)
    fused_aromatic_rings = set()
    for i, ring1 in enumerate(aromatic_rings):
        for ring2 in aromatic_rings[i+1:]:
            # Check if rings share at least one bond (a sign of being fused)
            shared_atoms = set(ring1).intersection(ring2)
            if len(shared_atoms) > 0:
                fused_aromatic_rings.update(shared_atoms)
    
    if len(fused_aromatic_rings) == 0:
        return False, "No fused aromatic rings, not a polycyclic arene"
    
    return True, "Structure is consistent with a polycyclic aromatic hydrocarbon"