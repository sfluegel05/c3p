"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: iridoid monoterpenoid
"""

from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    An iridoid monoterpenoid typically consists of a cyclopentane ring fused to a six-membered
    oxygen-containing heterocycle (pyran ring), forming a bicyclic system.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Check if molecule has at least two rings
    if len(atom_rings) < 2:
        return False, "Molecule does not have at least two rings"

    # Find fused ring systems (rings that share at least two atoms)
    fused_rings = []
    for i in range(len(atom_rings)):
        for j in range(i+1, len(atom_rings)):
            ring_i = set(atom_rings[i])
            ring_j = set(atom_rings[j])
            shared_atoms = ring_i & ring_j
            if len(shared_atoms) >= 2:
                fused_rings.append((atom_rings[i], atom_rings[j]))

    if not fused_rings:
        return False, "No fused ring systems found"

    found_iridoid = False
    for ring1, ring2 in fused_rings:
        size1 = len(ring1)
        size2 = len(ring2)
        # Check for cyclopentane (5-membered ring) and pyran (6-membered ring with oxygen)
        if (size1 == 5 and size2 ==6) or (size1 ==6 and size2 ==5):
            ring1_atoms = [mol.GetAtomWithIdx(idx) for idx in ring1]
            ring2_atoms = [mol.GetAtomWithIdx(idx) for idx in ring2]
            # Check if one ring is oxygen-containing heterocycle
            ring1_heteroatoms = [atom for atom in ring1_atoms if atom.GetAtomicNum() != 6]
            ring2_heteroatoms = [atom for atom in ring2_atoms if atom.GetAtomicNum() != 6]
            if (len(ring1_heteroatoms) == 1 and ring1_heteroatoms[0].GetAtomicNum() == 8) or \
               (len(ring2_heteroatoms) == 1 and ring2_heteroatoms[0].GetAtomicNum() == 8):
                found_iridoid = True
                break

    if not found_iridoid:
        return False, "No cyclopentane fused to oxygen-containing six-membered ring found"

    # Verify monoterpenoid backbone (10 carbon atoms)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 10:
        return False, f"Too few carbons for monoterpenoid (found {num_carbons}, need at least 10)"

    return True, "Contains fused cyclopentane and oxygen-containing six-membered ring typical of iridoid monoterpenoids"