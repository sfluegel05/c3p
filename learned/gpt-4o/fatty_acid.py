"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Count carbon atoms to determine chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4 or c_count > 28:
        return False, f"Carbon chain length of {c_count}, must be between 4 and 28"

    # Check if the molecule has aromatic rings
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings, not a fatty acid"

    # Check if all carbons are part of a non-cyclic, mostly linear chain (aliphatic)
    # The GetSSSR function returns all the smallest sets of smallest rings in the molecule.
    if Chem.GetSSSR(mol) > 0:
        return False, "Contains cyclic structures, not purely aliphatic"

    # If all checks pass, classify as fatty acid
    return True, "Contains aliphatic carbon chain with a carboxylic acid group"

# Example usage
# smiles = "CCCCC(=O)O"  # A simple example of a fatty acid, pentanoic acid
# print(is_fatty_acid(smiles))