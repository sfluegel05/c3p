"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: azole
Defined as any monocyclic heteroarene containing a five‐membered aromatic ring
with nitrogen (and possibly other non‐carbon atoms such as S or O) in the ring.
"""
from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is defined here as a compound containing a five-membered aromatic ring
    with at least one nitrogen atom, where the ring atoms are chosen among C, N, O and S.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains an azole ring, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "No rings found in the molecule"
    
    # Allowed atomic numbers in an azole ring: Carbon (6), Nitrogen (7), Oxygen (8), Sulfur (16)
    allowed_atoms = {6, 7, 8, 16}
    
    # Iterate over each ring (as a tuple of atom indices)
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue  # We need five-membered rings
        
        # Check that all atoms in the ring are aromatic and allowed
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            # Skip rings that are not fully aromatic
            continue
        
        # Check that all atoms in the ring are in the allowed set
        if not all(atom.GetAtomicNum() in allowed_atoms for atom in ring_atoms):
            continue
        
        # Check if at least one atom is nitrogen
        if not any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            continue
        
        # If we have reached this point, the ring qualifies as an azole ring
        return True, "Found a five-membered aromatic heterocycle containing nitrogen (azole ring)"
    
    # If none of the rings qualify
    return False, "No qualifying five-membered aromatic heterocycle (azole) found in the molecule"

# Example usage (comment out or remove in production)
if __name__ == "__main__":
    # Test with two example SMILES strings: one with an azole and one without.
    test_smiles_azole = "O1C(=NC=C1)CCCCC"  # 2-pentyloxazole from the examples
    test_smiles_non_azole = "CC(C)CC"         # A simple alkane
    
    result, reason = is_azole(test_smiles_azole)
    print("Test azole:", result, reason)
    
    result, reason = is_azole(test_smiles_non_azole)
    print("Test non-azole:", result, reason)