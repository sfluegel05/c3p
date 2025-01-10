"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene consists of multiple fused aromatic rings.
    
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
    
    # Detect aromatic rings using SMARTS for benzene-like rings
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')
    aromatic_benzenes = mol.GetSubstructMatches(benzene_pattern)
    
    # Check if there are at least two aromatic benzene rings
    if len(aromatic_benzenes) < 2:
        return False, f"Found {len(aromatic_benzenes)} aromatic benzene rings, need at least 2"
    
    # Check if rings are fused by deciding if they share atoms
    ring_info = mol.GetRingInfo()
    num_fused_rings = 0

    for i, ring1 in enumerate(ring_info.AtomRings()):
        for ring2 in ring_info.AtomRings():
            if ring1 != ring2 and any(atom in ring1 for atom in ring2):
                num_fused_rings += 1
                break

    if num_fused_rings < 2:
        return False, "Rings need to be fused"

    # Ensure no heteroatoms in rings (only carbon/hydrogen allowed for simple PAHs)
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() not in [6, 1] for idx in ring):
            return False, "Contains non-carbon, non-hydrogen atoms in aromatic rings"

    return True, "Molecule is a polycyclic arene with multiple aromatic rings fused together"

# Example usage
# print(is_polycyclic_arene("c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c21"))  # Should return True with reason