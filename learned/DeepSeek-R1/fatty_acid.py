"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is an aliphatic monocarboxylic acid with a chain of 4-28 carbons.

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

    # Check for exactly one carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups, need exactly 1"

    # Check for aromatic rings
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"

    # Check carbon chain length (total carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4 or c_count > 28:
        return False, f"Carbon count {c_count} not in 4-28 range"

    return True, "Aliphatic monocarboxylic acid with 4-28 carbons"