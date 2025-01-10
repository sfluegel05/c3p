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

    # Identify the carboxylic acid substructure
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"

    # To measure and verify the longest continuous carbon chain
    def is_linear_and_saturated(atom_idx, seen):
        """
        Verifies the longest continuous chain ensuring linear (no branch) and saturated (all single bonds).
        """
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6 or atom_idx in seen:
            return 0
        seen.add(atom_idx)
        chain_length = 1
        carbon_neighbors = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors += 1
                if atom.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx()).GetBondType() != Chem.rdchem.BondType.SINGLE:
                    return 0  # Unsaturated bond detected
            if carbon_neighbors > 2:
                return 0  # Branching detected
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in seen:
                chain_length += is_linear_and_saturated(neighbor.GetIdx(), seen)
        return chain_length

    # Start traversal from carboxylic acid group and visit each connected carbon atom
    max_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Only consider carbon atoms
            seen = set()
            chain_length = is_linear_and_saturated(atom.GetIdx(), seen)
            max_chain_length = max(max_chain_length, chain_length)
    
    # Check the chain length validity 
    if max_chain_length < 4 or max_chain_length > 36:
        return False, f"{max_chain_length} carbons found; expected a chain length between 4 and 36"
    
    return True, f"Contains {max_chain_length} carbons in a straight-chain saturated fatty acid structure"