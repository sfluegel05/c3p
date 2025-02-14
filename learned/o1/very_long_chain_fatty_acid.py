"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: very long-chain fatty acid
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is a fatty acid which has a chain length greater than or equal to C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid groups (C(=O)O[H])
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H1]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    # Identify ester functional groups (excluding carboxylic acid)
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H0][CX4]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if ester_matches:
        return False, "Ester functional group(s) found"

    # Identify amide functional groups
    amide_pattern = Chem.MolFromSmarts('[CX3](=O)[NX3][CX3]')
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if amide_matches:
        return False, "Amide functional group(s) found"

    # Count total number of carbons in molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if c_count >= 22:
        return True, f"Total carbon count is {c_count}, which is greater than or equal to 22"
    else:
        return False, f"Total carbon count is {c_count}, which is less than 22"