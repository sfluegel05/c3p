"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is typically a carboxylic acid with an aliphatic chain of fewer than 6 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group pattern - more general pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid functional group found"
    
    # Count the number of carbon atoms in the aliphatic chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Remove one for the carboxylic carbon
    aliphatic_chain_length = carbon_count - 1
    
    # SCFA should have an aliphatic chain of 5 or fewer carbons
    if aliphatic_chain_length > 5:
        return False, f"Aliphatic chain too long ({aliphatic_chain_length}), must be 5 or fewer excluding carboxyl carbon"

    return True, "Contains carboxylic acid group with a short aliphatic chain"