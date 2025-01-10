"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import deque

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester contains a lauric acid (C12:0) component attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    def count_carbon_chain(start_atom_idx):
        """Helper function to count longest carbon chain from starting atom"""
        visited = set()
        queue = deque([(start_atom_idx, 1)])  # (atom_idx, chain_length)
        max_length = 1
        
        while queue:
            current_idx, length = queue.popleft()
            current_atom = mol.GetAtomWithIdx(current_idx)
            visited.add(current_idx)
            
            # Look at neighbors
            for neighbor in current_atom.GetNeighbors():
                if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:  # Carbon
                    queue.append((neighbor.GetIdx(), length + 1))
                    max_length = max(max_length, length + 1)
        
        return max_length

    # For each ester group found
    for match in ester_matches:
        carbonyl_carbon = match[1]  # The C(=O) carbon
        
        # Count carbons in the chain starting from carbonyl carbon
        chain_length = count_carbon_chain(carbonyl_carbon)
        
        # Check if we found a 12-carbon chain
        if chain_length == 12:
            return True, "Contains dodecanoate ester group with 12-carbon chain"
            
    return False, "No dodecanoate ester group with appropriate chain length found"