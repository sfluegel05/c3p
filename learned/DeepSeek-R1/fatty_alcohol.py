"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol with a chain of 3 to over 27 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Find hydroxyl groups on non-aromatic, non-ring carbons
    hydroxyl_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1 and atom.GetTotalNumHs() >= 1:
            neighbor = atom.GetNeighbors()[0]
            if (neighbor.GetAtomicNum() == 6 and 
                not neighbor.GetIsAromatic() and 
                not neighbor.IsInRing()):
                hydroxyl_carbons.append(neighbor)
    
    if not hydroxyl_carbons:
        return False, "No hydroxyl group on aliphatic non-ring carbon"
    
    # Function to calculate longest aliphatic chain from a starting atom
    def get_chain_length(start_atom):
        visited = set()
        max_length = 1
        stack = [(start_atom, 1)]
        
        while stack:
            atom, length = stack.pop()
            visited.add(atom.GetIdx())
            
            for neighbor in atom.GetNeighbors():
                if (neighbor.GetAtomicNum() == 6 and
                    not neighbor.GetIsAromatic() and
                    not neighbor.IsInRing() and
                    neighbor.GetIdx() not in visited):
                    
                    new_length = length + 1
                    if new_length > max_length:
                        max_length = new_length
                    stack.append((neighbor, new_length))
        
        return max_length
    
    # Check chain lengths for all hydroxyl-bearing carbons
    max_chain = max(get_chain_length(carbon) for carbon in hydroxyl_carbons)
    
    if max_chain < 3:
        return False, f"Longest hydroxyl chain: {max_chain} (<3)"
    if max_chain > 27:
        return True, f"Aliphatic alcohol with {max_chain}-carbon chain (over 27)"
    
    return True, f"Aliphatic alcohol with {max_chain}-carbon chain"