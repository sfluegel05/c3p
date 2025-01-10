"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    This includes having one or more branching points in a hydrocarbon chain with a terminal carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    
    # Check for terminal carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Attempt to identify branched hydrocarbon chain
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing(): # Only consider non-ring carbons for branching
            carbon_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            # More than 2 carbon neighbors indicates potential branching
            if carbon_neighbors > 2:
                return True, "Contains a branched carbon chain with a terminal carboxylic acid group"
    
    # If no branching found, double check if could be due to stereochemistry or known complex structures (options beyond basic checks)
    return False, "No suitable branching pattern found in the straight chain structure leading to a carboxylic group"