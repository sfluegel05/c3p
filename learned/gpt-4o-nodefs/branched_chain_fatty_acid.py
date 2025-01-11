"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    Branched-chain fatty acids contain one or more carbon branching points.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for branching: atoms with more than two carbon atoms attached (non-linear)
    branched = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbors = atom.GetNeighbors()
            # Count neighboring carbon atoms
            carbon_neighbors = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 6)
            if carbon_neighbors > 2:  # More than 2 indicates branching
                branched = True
                break
    
    if not branched:
        return False, "No branching found in carbon chain"

    return True, "Contains carboxylic acid group and branch points in hydrocarbon chain"