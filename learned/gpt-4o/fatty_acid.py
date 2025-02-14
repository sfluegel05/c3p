"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for carbon atoms linked mainly in a linear chain
    c_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Count Carbon atoms
            # Check atom not in a ring, ensuring it's part of the main chain
            if not atom.IsInRing():
                c_count += 1

    # Natural fatty acids usually range from C4 to C28
    if not (4 <= c_count <= 28):
        return False, f"Carbon chain length of {c_count}, must be between 4 and 28"

    # Ensure molecule is mostly aliphatic (no aromatic rings)
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings, not a fatty acid"

    # Verify it is non-cyclic with no aromatic contents
    if Chem.rdmolops.GetSSSR(mol) > 0:
        return False, "Contains cyclic structures, not purely aliphatic"

    # Potentially allow additional functionalities or modifications (e.g. epoxides, hydroxy)
    # We may recognize specific features to support more complex fatty acids
    
    # Passed all checks, classifies as a fatty acid
    return True, "Contains aliphatic carbon chain with a carboxylic acid group"

# Example usage for testing
# smiles_ex = "CCCCC(=O)O"  # Pentanoic acid
# print(is_fatty_acid(smiles_ex))