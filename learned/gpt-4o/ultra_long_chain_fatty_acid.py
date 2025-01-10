"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    Defined as a fatty acid with a chain length greater than 27 carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find the carboxylic acid carbon (COOH group)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not matches:
        return False, "No carboxyl group found"

    carboxyl_carbon = matches[0][0]

    # Function to calculate the length of carbon chain starting from a given atom
    def get_chain_length_from(start_atom_idx, visited):
        to_visit = [start_atom_idx]
        length = 0
        
        while to_visit:
            atom_idx = to_visit.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                length += 1
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor.GetSymbol() == 'C':
                        to_visit.append(neighbor_idx)
        
        return length
    
    # Explore all chains starting from each carbon attached to the carboxyl group
    max_chain_length = 0
    visited = set()
    carboxyl_neighbors = mol.GetAtomWithIdx(carboxyl_carbon).GetNeighbors()
    for neighbor in carboxyl_neighbors:
        if neighbor.GetSymbol() == 'C':
            chain_length = get_chain_length_from(neighbor.GetIdx(), visited.copy())
            max_chain_length = max(max_chain_length, chain_length)
    
    if max_chain_length > 27:
        return True, f"Contains {max_chain_length} carbon atoms in chain, qualifies as ultra-long-chain fatty acid"
    else:
        return False, f"Contains {max_chain_length} carbon atoms in chain, does not qualify as ultra-long-chain fatty acid"