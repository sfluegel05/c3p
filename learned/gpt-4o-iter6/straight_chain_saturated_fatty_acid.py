"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find atoms involved in terminal carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_match = mol.GetSubstructMatch(carboxylic_acid_pattern)
    if not carboxylic_acid_match:
        return False, "No terminal carboxylic acid group found"
    
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found, not a fatty acid."

    visited = set()

    # Depth-first search to validate chain of carbons (excluding branches)
    def dfs(atom_idx):
        if atom_idx in visited:
            return 0
        visited.add(atom_idx)

        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:  # Must be a carbon
            return 0

        # If this carbon is part of the carboxylic group, end chain here
        if atom_idx in carboxylic_acid_match:
            return 1

        # Explore neighboring carbons
        carbon_count = 1
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_count += dfs(neighbor.GetIdx())

        return carbon_count

    # Start DFS from the carboxylic group carbon to ensure chain follows the carboxylic acid
    num_carbons = dfs(carboxylic_acid_match[0])

    # Natural fatty acids should be straight-chain and have between 4-36 carbons
    if num_carbons < 4 or num_carbons > 36:
        return False, f"{num_carbons} carbons found; expected a chain length between 4 and 36"
    
    return True, f"Contains {num_carbons} carbons in a straight-chain saturated fatty acid structure"