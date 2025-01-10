"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.

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

    # Look for carboxylic acid group pattern (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid functional group found"
    
    # Count the number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # SCFA should have chain of less than 6 carbons (including carboxyl carbon)
    if carbon_count > 6:
        return False, f"Too many carbon atoms ({carbon_count}), must be 6 or fewer"

    return True, "Contains carboxylic acid group with a short aliphatic chain"