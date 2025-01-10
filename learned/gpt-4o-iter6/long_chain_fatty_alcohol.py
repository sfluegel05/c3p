"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a carbon chain length ranging from C13 to C22 with
    at least one hydroxyl (-OH) group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for hydroxyl group (-OH) presence
    oh_group = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(oh_group):
        return False, "No hydroxyl group found"
    
    # Function to find the longest chain of carbon atoms using depth-first search
    def longest_carbon_chain(atom, visited, chain_length=0):
        atom_id = atom.GetIdx()
        visited.add(atom_id)
        max_length = chain_length
        
        # Iterate over neighbors
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:  # Only consider carbon atoms
                max_length = max(max_length, longest_carbon_chain(neighbor, visited, chain_length + 1))
        
        visited.remove(atom_id)
        return max_length
    
    # Find the maximum carbon chain length
    max_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Check only carbon atoms
            visited = set()
            max_chain_length = max(max_chain_length, longest_carbon_chain(atom, visited))
    
    # Verify if the longest carbon chain meets the C13 to C22 requirement
    if 13 <= max_chain_length <= 22:
        return True, f"Contains a carbon chain of length {max_chain_length} and a hydroxyl group"
    else:
        return False, f"Carbon chain length of {max_chain_length}, required between 13 and 22"

# Test the function with an example
print(is_long_chain_fatty_alcohol("CCCCCCCCCCCCCCCO"))  # Expected: True, reason: appropriate chain length and hydroxyl group