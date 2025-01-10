"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: Long-chain fatty alcohol (C13 to C22)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol is defined as a fatty alcohol with a chain length ranging from C13 to C22.

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

    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"

    # Function to calculate the longest carbon chain
    def get_longest_carbon_chain(mol):
        chains = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Carbon atom
                visited = set()
                stack = [(atom, 1)]  # (atom, chain_length)
                while stack:
                    current_atom, length = stack.pop()
                    visited.add(current_atom.GetIdx())
                    chains.append(length)
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                            stack.append((neighbor, length + 1))
        return max(chains) if chains else 0

    # Calculate the longest carbon chain
    longest_chain = get_longest_carbon_chain(mol)
    if longest_chain < 13 or longest_chain > 22:
        return False, f"Chain length is {longest_chain}, must be between 13 and 22"

    # Check if the hydroxyl group is attached to the longest chain
    hydroxyl_attached_to_longest_chain = False
    for match in hydroxyl_matches:
        hydroxyl_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in hydroxyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                # Check if this carbon is part of the longest chain
                visited = set()
                stack = [(neighbor, 1)]
                while stack:
                    current_atom, length = stack.pop()
                    visited.add(current_atom.GetIdx())
                    if length == longest_chain:
                        hydroxyl_attached_to_longest_chain = True
                        break
                    for next_neighbor in current_atom.GetNeighbors():
                        if next_neighbor.GetAtomicNum() == 6 and next_neighbor.GetIdx() not in visited:
                            stack.append((next_neighbor, length + 1))
                if hydroxyl_attached_to_longest_chain:
                    break
        if hydroxyl_attached_to_longest_chain:
            break

    if not hydroxyl_attached_to_longest_chain:
        return False, "Hydroxyl group not attached to the longest carbon chain"

    # Check for other functional groups that disqualify it as a fatty alcohol
    disallowed_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),  # Carboxylic acid
        Chem.MolFromSmarts("[CX3](=O)[OX2H0]"),  # Ester
        Chem.MolFromSmarts("[NX3]"),             # Amines
        Chem.MolFromSmarts("[SX2]"),             # Sulfides
    ]
    
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains disallowed functional groups"

    return True, f"Contains a hydroxyl group attached to a carbon chain of length {longest_chain}"