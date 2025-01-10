"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA consists of a coenzyme A moiety, a thioester linkage, and a short-chain acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more comprehensive CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)CCNC(=O)[C@H](O)COP(=O)(O)O[C@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No complete Coenzyme A moiety found"

    # Look for thioester group pattern, C(=O)SC with a fatty acyl chain of length 2 to 5 carbons
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)

    for match in thioester_matches:
        thioester_carbon = match[0]  # Gets the index of the ester carbon
        # Traverse neighbors to find a valid carbon chain length as per short-chain criteria
        carbon_chain_length = 0
        visited = set()
        queue = [thioester_carbon]
        
        while queue:
            atom_idx = queue.pop(0)
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom
                carbon_chain_length += 1
                if carbon_chain_length > 5:  # Check if the chain is still short
                    break
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:  # Ensure neighbor is carbon
                        queue.append(neighbor_idx)
        
        if 2 <= carbon_chain_length <= 5:
            return True, "Contains CoA moiety and short-chain fatty acyl thioester linkage"

    return False, "Does not satisfy short-chain fatty acyl criteria, or invalid linkages"