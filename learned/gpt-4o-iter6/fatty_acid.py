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
    if c_count < 4 or c_count > 28:
        return False, f"Carbon chain length {c_count} not in [4, 28]"

    # Check for presence of rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring structure(s), usually not characteristic of typical fatty acids"

    return True, "Valid fatty acid: Aliphatic monocarboxylic acid with between 4 and 28 carbons"