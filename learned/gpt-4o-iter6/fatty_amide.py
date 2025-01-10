"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is defined by an amide group derived from a long-chain fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of amide bond (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"
    
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Check for long carbon chains linked to amide bond
    for match in amide_matches:
        carbonyl_c_idx = match[0]  # Carbonyl Carbon
        nitrogen_idx = match[1]    # Nitrogen
        
        visited = set()
        carbon_chain_length = 0
        
        # Helper function for BFS search
        def bfs_carbon_chain(start_atom):
            queue = [(start_atom, 0)]  # (atom, depth)
            visited.add(start_atom.GetIdx())
            c_count = 0
            
            while queue:
                current_atom, depth = queue.pop(0)
                if current_atom.GetAtomicNum() == 6:
                    c_count += 1
                
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                        queue.append((neighbor, depth + 1))
                        visited.add(neighbor.GetIdx())

            return c_count
        
        # Count chain length from two sides (carbonyl and nitrogen)
        carbon_chain_length += bfs_carbon_chain(mol.GetAtomWithIdx(carbonyl_c_idx))
        carbon_chain_length += bfs_carbon_chain(mol.GetAtomWithIdx(nitrogen_idx))
        
        # Verify if enough carbon atoms are found
        if carbon_chain_length >= 12:  # Consider the threshold should capture typical long fatty acid chains
            return True, "Contains amide bond derived from a long chain fatty acid"

    return False, "Insufficient carbon chain length beyond amide bond"