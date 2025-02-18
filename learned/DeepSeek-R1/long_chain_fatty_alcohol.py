"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: long-chain fatty alcohol (C13-C22)
"""
from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol (C13-C22 chain).
    Requirements:
    1. Contains at least one hydroxyl group attached to an aliphatic carbon
    2. Hydroxyl-bearing carbon is part of a chain of 13-22 carbons connected by single bonds
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find all hydroxyl-bearing aliphatic carbons
    oh_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            if atom.GetTotalNumHs() == 1 and atom.GetDegree() == 1:  # -OH group
                neighbor = atom.GetNeighbors()[0]
                if neighbor.GetAtomicNum() == 6 and not neighbor.GetIsAromatic():
                    oh_carbons.append(neighbor)

    if not oh_carbons:
        return False, "No hydroxyl group on aliphatic carbon"

    # Calculate maximum chain length for each OH-bearing carbon
    def get_max_chain_length(start_atom):
        visited = set()
        max_length = 1
        stack = [(start_atom, 1)]
        
        while stack:
            atom, length = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            
            if length > max_length:
                max_length = length
                
            # Explore all adjacent aliphatic carbons connected by single bonds
            for neighbor in atom.GetNeighbors():
                if (neighbor.GetAtomicNum() == 6 and 
                    not neighbor.GetIsAromatic() and 
                    atom.GetBondTo(neighbor).GetBondType() == Chem.BondType.SINGLE):
                    stack.append((neighbor, length + 1))
        
        return max_length

    # Check chain lengths for all hydroxyl-bearing carbons
    for carbon in oh_carbons:
        chain_length = get_max_chain_length(carbon)
        if 13 <= chain_length <= 22:
            return True, f"Hydroxyl on C{chain_length} aliphatic chain"
    
    return False, "No hydroxyl on C13-C22 chain"