"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a straight or branched carbon chain ranging from C13 to C22 with
    at least one terminal hydroxyl (-OH) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxyl group (-OH) presence
    oh_group = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(oh_group):
        return False, "No hydroxyl group found"

    # Calculate carbon chain lengths
    chain_lengths = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Only consider carbon atoms
            visited = set()
            visited.add(atom.GetIdx())
            chain_lengths.append(longest_carbon_chain_including_oh(atom, visited, mol))

    # Verify if the longest carbon chain meets the C13 to C22 requirement
    if not chain_lengths:
        return False, "No valid carbon chain found"
        
    max_chain_length = max(chain_lengths)
    if 13 <= max_chain_length <= 22:
        return True, f"Contains a carbon chain of length {max_chain_length} and a hydroxyl group"
    else:
        return False, f"Carbon chain length of {max_chain_length}, required between 13 and 22"

def longest_carbon_chain_including_oh(atom, visited, mol, chain_length=1):
    max_length = chain_length
    for bond in atom.GetBonds():
        neighbor = bond.GetOtherAtom(atom)
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:  # Continue with carbon
            visited.add(neighbor.GetIdx())
            max_length = max(max_length, longest_carbon_chain_including_oh(neighbor, visited, mol, chain_length + 1))
            visited.remove(neighbor.GetIdx())
    return max_length

# Test the function with an example
print(is_long_chain_fatty_alcohol("CCCCCCCCCCCCCCCO"))  # Expected: True, appropriate chain length and hydroxyl group