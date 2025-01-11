"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: Long-chain fatty acid (C13 to C22)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid is defined as having a carbon chain length between 13 and 22 carbons,
    with at least one carboxylic acid group (-COOH). The molecule may contain other functional groups
    and cyclic structures as long as the main carbon chain length is within the specified range.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13 or c_count > 22:
        return False, f"Carbon chain length {c_count} is not between 13 and 22"

    # Check if the molecule has a reasonable number of rotatable bonds for a fatty acid
    # This helps exclude very rigid structures that are unlikely to be fatty acids
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for a long-chain fatty acid"

    # If all checks pass, classify as a long-chain fatty acid
    return True, f"Contains a carboxylic acid group and a carbon chain length of {c_count} (C13 to C22)"