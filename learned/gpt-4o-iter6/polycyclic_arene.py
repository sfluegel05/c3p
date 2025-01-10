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
    
    # Get info about rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = []
    
    # Find all aromatic rings
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(set(ring))
    
    # Check if there are at least two aromatic rings
    if len(aromatic_rings) < 2:
        return False, f"Molecule has only {len(aromatic_rings)} aromatic rings, insufficient for polycyclic arene"
    
    # Check for fused aromatic rings
    fused_ring_count = 0
    for i, ring1 in enumerate(aromatic_rings):
        for ring2 in aromatic_rings[i+1:]:
            # If there's at least one shared atom, the rings are fused
            if ring1.intersection(ring2):
                fused_ring_count += 1
                break
    
    # Determine if it is a polycyclic aromatic system
    if fused_ring_count >= 1:
        # Ensure it primarily consists of carbon (as a PAH)
        if all(atom.GetAtomicNum() in [6] for atom in mol.GetAtoms() if atom.GetIsAromatic()):
            return True, "Molecule is a polycyclic arene with multiple fused aromatic rings"
    
    return False, f"Molecule does not meet polycyclic arene criteria with {len(aromatic_rings)} aromatic rings and {fused_ring_count} fused rings"