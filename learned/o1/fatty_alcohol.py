"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
    Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove hydrogens for efficiency
    mol = Chem.RemoveHs(mol)
    
    # Identify all carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]

    if not carbon_atoms:
        return False, "No carbon atoms found"

    # Find all chains starting from carbon atoms
    from collections import deque

    longest_chain_length = 0
    longest_chain = []

    for atom in carbon_atoms:
        visited = set()
        queue = deque()
        queue.append([atom.GetIdx()])
        
        while queue:
            path = queue.popleft()
            current_atom_idx = path[-1]
            current_atom = mol.GetAtomWithIdx(current_atom_idx)
            visited.add(current_atom_idx)
            
            # Extend path
            neighbors = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited]
            for nbr in neighbors:
                new_path = path + [nbr.GetIdx()]
                queue.append(new_path)
                if len(new_path) > longest_chain_length:
                    longest_chain_length = len(new_path)
                    longest_chain = new_path

    if longest_chain_length < 3:
        return False, f"Longest carbon chain has {longest_chain_length} carbons, less than 3"

    # Check for hydroxyl group attached to any carbon in the longest chain
    found = False
    for idx in longest_chain:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 1:
                    # Check that oxygen is hydroxyl group (OH)
                    if nbr.GetDegree() == 1 and not nbr.IsInRing():
                        found = True
                        break
            if found:
                break

    if not found:
        return False, "No hydroxyl group attached to carbon in the longest chain"

    return True, f"Contains aliphatic alcohol with a chain of {longest_chain_length} carbons"

# Example usage:
# smiles = "CCCCCCCCCCCCCC(O)CCC"  # tetradecan-5-ol
# result, reason = is_fatty_alcohol(smiles)
# print(result, reason)