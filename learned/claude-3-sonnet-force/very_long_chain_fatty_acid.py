"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:36976 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX1-,N]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carboxylic acid group found"

    # Get the longest linear chain length
    chain_length = Chem.Lipinski.GetLongestChain(mol)

    # Check chain length
    if chain_length <= 22:
        return False, f"Chain length {chain_length} is too short for very long-chain fatty acid"

    # Check for linear chain (no rings, no substituents other than -C(=O)O)
    linear_pattern = Chem.MolFromSmarts("[C][C](=[O])[O]")
    matches = mol.GetSubstructMatches(linear_pattern)
    if len(matches) != 1:
        return False, "Not a linear chain or has additional substituents"

    # Check for additional rings
    if Chem.MolFromSmiles(smiles).GetRingInfo().NumRings() > 0:
        return False, "Contains rings, not a pure fatty acid"

    # Check for additional functional groups
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    if mol.HasSubstructMatch(alcohol_pattern):
        return False, "Contains alcohol groups, not a pure fatty acid"

    return True, f"Linear fatty acid with chain length {chain_length} (> C22)"