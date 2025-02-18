"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A specific fragments: check for the entire CoA scaffold
    coA_pattern = Chem.MolFromSmarts("OC(=O)CNC(=O)[C@H](O)[C@H](C)C(=O)NCCC(=O)NCC")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Coenzyme A structure not found"

    # Thioester linkage: S-C(=O)
    thioester_pattern = Chem.MolFromSmarts("S-[CX3](=O)-")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Check if the chain length starting from the carbonyl carbon is 2 to 6 carbons
    for match in thioester_matches:
        carbonyl_carbon_idx = match[1]

        # Start BFS from the carbon atom next in the chain
        neighbors = [a for a in mol.GetAtomWithIdx(carbonyl_carbon_idx).GetNeighbors() if a.GetIdx() != match[0]]
        
        chain_length = 0
        queue = [(n, 1) for n in neighbors]  # Store (atom, current_chain_length)
        visited = set()

        while queue:
            current_atom, current_length = queue.pop(0)
            current_idx = current_atom.GetIdx()

            if current_idx in visited:
                continue
            
            visited.add(current_idx)

            if current_atom.GetAtomicNum() == 6:
                chain_length = current_length

                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited:
                        queue.append((neighbor, current_length + 1))

        # Validate if chain length is 2 to 6
        if 2 <= chain_length <= 6:
            return True, "Contains Coenzyme A linked with a short-chain fatty acid"

    return False, "No short-chain fatty acid attached to CoA"