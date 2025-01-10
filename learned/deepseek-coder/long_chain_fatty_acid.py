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
    with a carboxylic acid group (-COOH).

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

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13 or c_count > 22:
        return False, f"Carbon chain length {c_count} is not between 13 and 22"

    # Check if the molecule is primarily a carbon chain with a carboxylic acid group
    # This is a heuristic and may not cover all edge cases
    # We look for a linear or branched carbon chain with the carboxylic acid group at one end
    # and no other significant functional groups
    # This is a simplified check and may need refinement for complex cases
    # For example, we exclude molecules with rings or multiple carboxylic acid groups
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Molecule contains rings, which are not typical for long-chain fatty acids"

    # Count the number of carboxylic acid groups
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) > 1:
        return False, "Multiple carboxylic acid groups found"

    # Check for other significant functional groups (e.g., alcohols, amines, etc.)
    # This is a heuristic and may need adjustment
    other_functional_groups = ["[OH]", "[NH2]", "[N]", "[S]", "[P]"]
    for group in other_functional_groups:
        pattern = Chem.MolFromSmarts(group)
        if mol.HasSubstructMatch(pattern):
            return False, f"Found unexpected functional group: {group}"

    # If all checks pass, classify as a long-chain fatty acid
    return True, f"Contains a carboxylic acid group and a carbon chain length of {c_count} (C13 to C22)"