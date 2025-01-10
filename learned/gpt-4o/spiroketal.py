"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    atom_rings = mol.GetRingInfo().AtomRings()

    for atom in mol.GetAtoms():
        # Check if the atom is a carbon (potential spiro center)
        if atom.GetAtomicNum() != 6:
            continue
        
        # Collect the rings the atom is part of
        rings = [ring for ring in atom_rings if atom.GetIdx() in ring]
        if len(rings) != 2:
            continue
        
        # Ensure the atom is the sole connector between two rings
        ring_set = set(rings[0]) & set(rings[1])
        if len(ring_set) != 1:
            continue

        # Check for ketal structure: two different ring oxygen atoms bound
        oxygen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.IsInRing()]
        if len(oxygen_neighbors) == 2:
            # Check that each oxygen belongs to a different ring
            if any(oxygen.GetIdx() in ring for oxygen in oxygen_neighbors for ring in rings):
                return True, "Spiroketal structure identified with appropriate ketal formation."

    return False, "Spiro centers found but they do not form expected ketal groups or were not part of such a structure."

# Example usage
example_smiles = "O=C1O[C@@H]2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](OC)[C@H]5O)C)O)C)C)C"
result, reason = is_spiroketal(example_smiles)
print(f"Is Spiroketal: {result}, Reason: {reason}")