"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is defined by an amide group derived from a long chain fatty acid.

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

    # Find the carbonyl group within the amide bond
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    for match in amide_matches:
        carbonyl_c_idx = match[0]
        
        # Perform a breadth-first search starting from the carbonyl carbon to ensure a long carbon chain
        visited = set()
        queue = [carbonyl_c_idx]
        carbon_chain_length = 0
        
        while queue:
            current_idx = queue.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            
            current_atom = mol.GetAtomWithIdx(current_idx)
            if current_atom.GetAtomicNum() == 6:  # Check if it's a carbon atom
                carbon_chain_length += 1
                if carbon_chain_length >= 8:  # Chain is long enough
                    return True, "Contains amide bond derived from a long chain fatty acid"

            # Add neighbors to the search queue
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    queue.append(neighbor_idx)

    return False, "Insufficient carbon chain length beyond amide bond"