"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene consists of multiple fused aromatic rings, typically hydrocarbon.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Verify if the molecule contains aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        # Check if all atoms in the ring are aromatic
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
    
    # We need at least two aromatic rings for polycyclic arenes
    if aromatic_ring_count < 2:
        return False, f"Molecule has only {aromatic_ring_count} aromatic rings, which is insufficient for polycyclic arene classification"
    
    # Check for fused rings based on shared bonds
    bond_rings = [set(ring) for ring in ring_info.BondRings()]
    fused_rings = 0
    for i, ring1 in enumerate(bond_rings):
        for ring2 in bond_rings[i+1:]:
            # If there is at least one shared bond, the rings are fused
            if ring1.intersection(ring2):
                fused_rings += 1
                break
    
    # Determine overall polycyclic aromatic system possibility
    if aromatic_ring_count >= 2 and fused_rings >= 1:
        return True, "Molecule is a polycyclic arene with multiple fused aromatic rings"
    
    return False, f"Molecule does not meet polycyclic arene criteria with {aromatic_ring_count} aromatic rings and {fused_rings} fused rings"

# Example usage
# print(is_polycyclic_arene("c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c21"))  # Should return True with reason