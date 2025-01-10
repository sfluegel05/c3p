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

    # Check for alcohol group (oxygen connected to a carbon, checking for both primary, secondary)
    alcohol_pattern_primary = Chem.MolFromSmarts('[CX4][OX2H]')
    alcohol_pattern_secondary = Chem.MolFromSmarts('[CX3][OX2H]')
    alcohol_pattern_tertiary = Chem.MolFromSmarts('[CX2][OX2H]')
    if not (mol.HasSubstructMatch(alcohol_pattern_primary) or 
            mol.HasSubstructMatch(alcohol_pattern_secondary) or 
            mol.HasSubstructMatch(alcohol_pattern_tertiary)):
        return False, "No proper alcohol group found"

    # Check for total carbon count (at least 14)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14:
        return False, "Not enough total carbons for a long hydrocarbon chain"

    # Success: Detected a long-chain fatty alcohol
    return True, "Contains long hydrocarbon chain with an alcohol group"

# Test the function with a sample SMILES string from the list:
print(is_long_chain_fatty_alcohol('CCCCCCCCCC/C=C/CCCCCCCCCCCCO'))  # Example: (13Z)-docosen-1-ol