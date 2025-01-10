"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing carboxylic acid group (C(=O)O)"
    
    # Check for long aliphatic chain, typically 8 to 28 carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if carbon_count < 8:
        return False, f"Too few carbons ({carbon_count}) for long-chain fatty acid"
    if carbon_count > 28:
        return False, f"Too many carbons ({carbon_count}), exceeds typical long-chain fatty acid length"

    return True, "Structure supports classification as a long-chain fatty acid"