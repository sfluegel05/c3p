"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem

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

    # Look for a hydroxyl group (-OH) attached to carbon
    hydroxyl_pattern = Chem.MolFromSmarts("CO")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group on a carbon chain found"

    # Function to count the longest continuous carbon chain
    def longest_chain_length(mol):
        max_length = 0
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                visited = set()
                current_chain = []
                current_chain.append(bond.GetBeginAtomIdx())
                current_chain.append(bond.GetEndAtomIdx())
                while current_chain:
                    atom_idx = current_chain.pop()
                    if atom_idx not in visited:
                        visited.add(atom_idx)
                        atom = mol.GetAtomWithIdx(atom_idx)
                        if atom.GetAtomicNum() == 6:  # Carbon
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                                    current_chain.append(neighbor.GetIdx())
                max_length = max(max_length, len(visited))
        return max_length

    longest_chain = longest_chain_length(mol)

    # Check if carbon count is within 13 to 22
    if longest_chain < 13 or longest_chain > 22:
        return False, f"Carbon chain length is {longest_chain}, expected between 13 and 22"
    
    return True, "Has a hydroxyl group and the carbon chain length is within range for long-chain fatty alcohol"