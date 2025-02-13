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

    # Identify carboxylic acid substructure
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"

    # To measure the longest straight chain
    def longest_straight_chain(atom_idx, seen):
        """
        Recursively finds the longest straight chain from a starting atom index.
        """
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6 or atom_idx in seen:
            return 0
        seen.add(atom_idx)
        chain_length = 1
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in seen:
                chain_length = max(chain_length, 1 + longest_straight_chain(neighbor.GetIdx(), seen))
        return chain_length
    
    # Consider all carbon atoms and their maximum connections
    max_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Only consider carbon atoms
            seen = set()
            chain_length = longest_straight_chain(atom.GetIdx(), seen)
            max_chain_length = max(max_chain_length, chain_length)
    
    # Check the chain length validity 
    if max_chain_length < 4 or max_chain_length > 36:
        return False, f"{max_chain_length} carbons found; expected a chain length between 4 and 36"
    
    return True, f"Contains {max_chain_length} carbons in a straight-chain saturated fatty acid structure"