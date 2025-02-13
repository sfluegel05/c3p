"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is defined as a molecule with a long aliphatic chain and a terminal carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the carboxylic acid group pattern (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for a long aliphatic carbon chain (saturated/unsaturated)
    # Long chain is defined as more than 5 carbon atoms
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCC")  # At least 6 contiguous carbon atoms
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long aliphatic chain (more than 5 Cs) found"

    return True, "Contains a carboxylic acid group and a long aliphatic chain"

# Testing the function with an example
example_smiles = "OC(=O)C/C=C\CCCC/C=C\CCCCCC"
result, reason = is_fatty_acid(example_smiles)
print(f"Is fatty acid: {result}, Reason: {reason}")