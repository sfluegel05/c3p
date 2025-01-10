"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: Long-chain fatty alcohol (C13 to C22)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol is defined as a fatty alcohol with a chain length ranging from C13 to C22.

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

    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"

    # Count the number of carbon atoms in the longest chain
    longest_chain = rdMolDescriptors.CalcLongestChain(mol)
    if longest_chain < 13 or longest_chain > 22:
        return False, f"Chain length is {longest_chain}, must be between 13 and 22"

    # Check for other functional groups that disqualify it as a fatty alcohol
    # For example, carboxylic acids, esters, etc.
    disallowed_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),  # Carboxylic acid
        Chem.MolFromSmarts("[CX3](=O)[OX2H0]"),  # Ester
        Chem.MolFromSmarts("[NX3]"),             # Amines
        Chem.MolFromSmarts("[SX2]"),             # Sulfides
    ]
    
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains disallowed functional groups"

    return True, f"Contains a hydroxyl group and a carbon chain of length {longest_chain}"