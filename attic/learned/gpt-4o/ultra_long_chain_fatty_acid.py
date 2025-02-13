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

    # Find carboxyl carbon (COOH group) and start counting from there
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not matches:
        return False, "No carboxyl group found"
    
    # Assuming the first match is the relevant carboxyl group
    carboxyl_carbon = matches[0][0]

    # Use BFS to count consecutive carbon atoms in the chain attached to the carboxyl carbon
    visited = set()
    to_visit = [carboxyl_carbon]
    carbon_count = 0
    
    while to_visit:
        atom_idx = to_visit.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != 'C':
            continue
        carbon_count += 1
        
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited:
                to_visit.append(neighbor_idx)

    if carbon_count > 27:
        return True, f"Contains {carbon_count} carbon atoms in chain, qualifies as ultra-long-chain fatty acid"
    else:
        return False, f"Contains {carbon_count} carbon atoms in chain, does not qualify as ultra-long-chain fatty acid"