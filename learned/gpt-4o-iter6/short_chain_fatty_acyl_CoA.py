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

    # Enhanced CoA moiety matching pattern (should account for key structural features of CoA)
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCSC(=O)CCNC(=O)[C@H](O)COP(=O)(O)O[C@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No complete Coenzyme A moiety found"

    # Thioester group with variable chain length, at least 2 to 5 carbons attached
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)

    for match in thioester_matches:
        thioester_carbon = match[0]  # Thioester carbon
        visited = set()
        carbon_queue = [(thioester_carbon, 0)]  # Store (atom index, chain length)

        while carbon_queue:
            atom_idx, chain_length = carbon_queue.pop(0)
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)

            if atom.GetAtomicNum() == 6:  # Carbon
                chain_length += 1
                if chain_length > 5:
                    break
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:  # Carbon
                        carbon_queue.append((neighbor_idx, chain_length))

        if 2 <= chain_length <= 5:  # Chain length must be short
            return True, "CoA moiety and short-chain fatty acyl thioester linkage present"

    return False, "Does not satisfy short-chain fatty acyl criteria, or invalid linkages"