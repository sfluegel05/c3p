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

    # Check for presence of hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group (-OH) found"

    # For each hydroxyl group, check for connected carbon chain of at least 3 carbons
    for match in hydroxyl_matches:
        o_idx = match[0]  # Index of the oxygen atom
        o_atom = mol.GetAtomWithIdx(o_idx)
        
        # Find the carbon atom attached to the hydroxyl oxygen
        c_neighbors = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not c_neighbors:
            continue  # No carbon neighbor found
        c_atom = c_neighbors[0]
        
        # Perform BFS to count connected carbons (including branches)
        visited = set()
        queue = [c_atom]
        carbon_count = 0
        while queue:
            atom = queue.pop(0)
            idx = atom.GetIdx()
            if idx in visited:
                continue
            visited.add(idx)
            carbon_count += 1
            # Traverse to neighboring carbon atoms
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    queue.append(neighbor)
        
        if carbon_count >= 3:
            return True, f"Molecule has a hydroxyl group connected to a carbon chain of {carbon_count} carbons"
    
    return False, "No hydroxyl group connected to a carbon chain of at least 3 carbons found"

__metadata__ = {
    'chemical_class': {
        'name': 'fatty alcohol',
        'definition': 'An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms. Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.',
    }
}