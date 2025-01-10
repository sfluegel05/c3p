"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    Long-chain fatty alcohols contain a long hydrocarbon chain (typically 14 or more carbons)
    with an alcohol (-OH) group at the end or within the chain.

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

    # Check for alcohol group (oxygen connected to a saturated carbon)
    alcohol_pattern = Chem.MolFromSmarts('[CX4H2,CX4H3][OX2H]')
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No proper alcohol group found"

    # Check for long carbon chain of at least 14 carbons
    carbon_chain_pattern = Chem.MolFromSmarts('[CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0]')
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long hydrocarbon chain found"

    return True, "Contains long hydrocarbon chain with an alcohol group"

# Test the function with a sample SMILES string from the list:
print(is_long_chain_fatty_alcohol('CCCCCCCCCC/C=C/CCCCCCCCCCCCO'))  # Example: (13Z)-docosen-1-ol