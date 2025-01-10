"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem

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
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Locate the carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Initialize count of carbons in aliphatic chains
    carbon_count = 0

    # Traverse atoms and count carbon atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            carbon_count += 1

    # Now determine if chain is longer than C22
    if carbon_count > 22:
        return True, f"Contains {carbon_count} carbons, qualifies as a very long-chain fatty acid"
    else:
        return False, f"Contains {carbon_count} carbons, does not exceed C22"

# Example usage:
# result, reason = is_very_long_chain_fatty_acid("C(O)(=O)CCCCCCCCCCCCCCCCCCCCCC(O)=O")
# print(result, reason)