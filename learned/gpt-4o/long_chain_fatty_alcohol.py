"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a primary alcohol group and a carbon chain length of 13 to 22 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for primary alcohol group (C-O-H) at the terminal position of carbon chain
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    if not mol.HasSubstructMatch(primary_alcohol_pattern):
        return False, "No primary alcohol group found"
    
    # Count the longest carbon chain
    carbon_chain_length = rdMolDescriptors.CalcMolLongestPath(mol)
    
    # Verify the carbon chain is between 13 and 22 carbons
    if carbon_chain_length < 13 or carbon_chain_length > 22:
        return False, f"Carbon chain length is {carbon_chain_length}, outside the 13-22 range"

    return True, f"Contains a primary alcohol group and has a carbon chain length of {carbon_chain_length}."