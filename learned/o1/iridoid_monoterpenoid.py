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
    oxygen-containing heterocycle.

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

    # Define SMARTS pattern for cyclopentane ring fused to a six-membered oxygen heterocycle
    # The pattern represents a bicyclic system with shared atoms between the rings
    fused_ring_pattern = Chem.MolFromSmarts("""
        [
            [#6]-1          # Carbon atom labeled 1
            ~[#6]~[#6]      # Connected to two other carbons
            ~[#6]-2         # Connected to a carbon labeled 2
            ~[#8]~[#6]~[#6] # Six-membered ring with an oxygen
            ~1~2            # Ring closures between atoms 1 and 2
        ]
    """)
    if fused_ring_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for the fused ring system
    if not mol.HasSubstructMatch(fused_ring_pattern):
        return False, "No cyclopentane fused to six-membered oxygen heterocycle found"

    # Verify monoterpenoid backbone (10 carbon atoms derived from isoprene units)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 10:
        return False, f"Too few carbons for monoterpenoid (found {num_carbons}, need at least 10)"

    return True, "Contains cyclopentane fused to six-membered oxygen heterocycle typical of iridoid monoterpenoids"