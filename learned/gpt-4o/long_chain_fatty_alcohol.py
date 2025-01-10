"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a chain length ranging from C13 to C22 with a hydroxyl group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a hydroxyl group (-OH) attached to aliphatic carbon
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")  # Aliphatic carbon with OH
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    if not hydroxyl_matches:
        return False, "No aliphatic hydroxyl group (-OH) found on a carbon chain"

    # Check longest aliphatic carbon chain with the attached hydroxyl group
    def max_aliphatic_chain_length(attached_idx):
        max_chain_length = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Carbon
                visited = set()
                stack = [(atom.GetIdx(), 1)]  # (Current atom index, chain length)
                while stack:
                    current_atom_idx, current_length = stack.pop()
                    if current_atom_idx not in visited:
                        visited.add(current_atom_idx)
                        current_atom = mol.GetAtomWithIdx(current_atom_idx)
                        for neighbor in current_atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                                stack.append((neighbor.GetIdx(), current_length + 1))
                        if current_atom.GetIdx() == attached_idx:
                            max_chain_length = max(max_chain_length, current_length)
        return max_chain_length

    # Evaluate each hydroxyl group's attachment
    valid_alcohol = False
    reason = ""
    for match in hydroxyl_matches:
        carbon_idx, oxygen_idx = match
        chain_length = max_aliphatic_chain_length(carbon_idx)
        if 13 <= chain_length <= 22:
            valid_alcohol = True
            reason = "Has a hydroxyl group and the carbon chain length is within range for long-chain fatty alcohol"
            break
        else:
            reason = f"Carbon chain length is {chain_length}, expected between 13 and 22"

    return valid_alcohol, reason