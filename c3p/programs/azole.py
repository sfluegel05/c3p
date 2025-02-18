"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: azole
Defined as any monocyclic heteroarene consisting of a five‐membered aromatic ring
containing at least one nitrogen atom. The ring must be composed only of C, N, O, and S,
and must not be fused to any other ring.
"""
from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule contains an azole ring (i.e. a monocyclic, isolated, five‐membered
    aromatic heterocycle with at least one nitrogen) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains an azole, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update aromaticity perception (if needed)
    Chem.SanitizeMol(mol)

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in the molecule"
    
    # Pre-calculate how many rings each atom is part of;
    # if an atom belongs to more than one ring, it means the ring is fused.
    atom_ring_count = {}
    for ring in atom_rings:
        for idx in ring:
            atom_ring_count[idx] = atom_ring_count.get(idx, 0) + 1

    # Allowed atomic numbers in an azole ring: C (6), N (7), O (8) and S (16)
    allowed_atoms = {6, 7, 8, 16}
    
    # Iterate over each ring
    for ring in atom_rings:
        # We only care about five-membered rings
        if len(ring) != 5:
            continue
        
        # Check if the ring is isolated (non-fused)
        # If any atom in the ring is part of >1 ring, skip this candidate ring.
        if any(atom_ring_count[idx] > 1 for idx in ring):
            continue
        
        # Get atoms in this ring
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Check that every atom in the ring is aromatic.
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        
        # Check that every atom in the ring is one of the allowed elements.
        if not all(atom.GetAtomicNum() in allowed_atoms for atom in ring_atoms):
            continue
        
        # Check if the ring contains at least one nitrogen atom.
        if not any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            continue
        
        # If we reach this point, the candidate ring qualifies as an azole ring.
        return True, "Found an isolated five-membered aromatic heterocycle containing nitrogen (azole ring)"
    
    # None of the rings qualified as an azole ring.
    return False, "No qualifying isolated five-membered aromatic heterocycle (azole) found in the molecule"


# Example usage (you can remove or comment out this test section in production)
if __name__ == "__main__":
    # Test with an azole: 2-Pentyloxazole as given in the examples.
    test_smiles_azole = "O1C(=NC=C1)CCCCC"
    result, reason = is_azole(test_smiles_azole)
    print("Test azole:", result, reason)
    
    # Test with a non-azole: simple alkane.
    test_smiles_non_azole = "CC(C)CC"
    result, reason = is_azole(test_smiles_non_azole)
    print("Test non-azole:", result, reason)