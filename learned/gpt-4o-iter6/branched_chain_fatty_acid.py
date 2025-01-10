"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid has a long carbon chain with one or more alkyl substituents and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (R-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for alkyl chain by counting carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # Assumes at least 10 carbons indicative of long chain
        return False, "Insufficient carbon count for alkyl chain"

    # Look for branching
    branch_pattern = Chem.MolFromSmarts("[CX4][CX4,CX3]([CX4,CX3])[CX4,CX3]")
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No branch points found"

    return True, "Contains sufficient carbon chain with branch points and a carboxylic acid group"