"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: fatty aldehyde (CHEBI:35581)
"""

from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde formally arising from reduction of the carboxylic acid group
    of its corresponding fatty acid, having a carbonyl group at one end of the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define aldehyde SMARTS pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H]=O")

    # Check for at least one aldehyde group
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if len(aldehyde_matches) == 0:
        return False, "No aldehyde group found."

    # Check each aldehyde group
    for match in aldehyde_matches:
        aldehyde_carbon_idx = match[0]
        aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_carbon_idx)
        
        # Get the non-oxygen neighbor atoms of the aldehyde carbon
        neighbors = [atom for atom in aldehyde_carbon.GetNeighbors() if atom.GetAtomicNum() != 8]
        
        # Check if the aldehyde carbon is connected to at least one carbon (part of a chain)
        chain_found = False
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:
                chain_found = True
                break
        if not chain_found:
            continue  # Check the next aldehyde group
        
        # Optional: traverse the carbon chain to check chain length (not strictly necessary)
        # For fatty aldehydes, chain length can vary, so we won't enforce a minimum
        
        return True, "Contains aldehyde group at the end of a carbon chain."

    # If no aldehyde groups are connected to a carbon chain
    return False, "Aldehyde group not connected to a carbon chain."