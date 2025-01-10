"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is characterized as an aliphatic monocarboxylic acid with a chain of 4 to 28 carbons 
    (usually unbranched and even-numbered), which may be saturated or unsaturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Get number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Fatty acids typically have between 4 and 28 carbon atoms
    # Adjust criteria to allow for more complexity than strictly chain length
    if c_count < 4:
        return False, f"Carbon chain length {c_count} too short for typical fatty acid"
    
    # We need to be cautious with long chains that might just be complex esters
    # Look at length but allow higher values with more checks
    if c_count > 28:
        larger_chain_info = f"Carbon chain length {c_count} for complexities beyond typical"
        # If structure is simple and exceeds typical complexity of these types
        if mol.GetRingInfo().NumRings() == 0:
            return False, larger_chain_info
        connected_heavy_atoms = mol.GetNumHeavyAtoms()
        if connected_heavy_atoms < c_count * 1.5: # Allow some connectivity above simple alkane chains
            return True, f"Larger but valid fatty acid: {larger_chain_info}"

    # Check for presence of rings and ensure they are small if present
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0 and all(len(r) > 6 for r in ring_info.AtomRings()):
        return False, "Contains large ring structure(s), not characteristic of typical fatty acids"

    return True, "Valid fatty acid: Aliphatic monocarboxylic acid with primarily aliphatic character"