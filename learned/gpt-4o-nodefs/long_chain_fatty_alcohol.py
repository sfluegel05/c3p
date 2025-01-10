"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    Long-chain fatty alcohols contain a long hydrocarbon chain (typically 14 or more carbons)
    with an alcohol (-OH) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for alcohol group (O connected directly to a carbon that is not part of carboxylic acid)
    alcohol_pattern = Chem.MolFromSmarts('[CX4][OX2H]')
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"
    
    # Check for long hydrocarbon chain of at least 14 carbons
    carbon_chain_pattern = Chem.MolFromSmarts("[C;R0](=O=[O+])([C;R0](*[OX2H]));!*[!#6]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)

    # Check if there is any chain with at least 14 carbon atoms
    long_chain_found = any(len(match) >= 14 for match in carbon_chain_matches)
    
    if not long_chain_found:
        return False, "No long hydrocarbon chain found"
    
    return True, "Contains long hydrocarbon chain with an alcohol group"