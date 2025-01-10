"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Attempt to identify the Coenzyme A pattern
    # CoA moiety often includes a core structure but can have flexible ends; using a broader pattern
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)COP(=O)(O)O[C@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No complete Coenzyme A moiety found"

    # Improve the detection of thioester group with short-chain acyl group
    # The thioester pattern targets the sulfur linkage to an acyl group
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)

    for match in thioester_matches:
        thioester_carbon = match[0]  # Thioester carbon
        visited = set()
        chain_length = 0
        atom_queue = [(thioester_carbon)]
        
        while atom_queue:
            atom_idx = atom_queue.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)

            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                chain_length += 1
                if chain_length > 5:  # Exceeded short-chain length
                    break
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited:  # Check all neighbors
                        atom_queue.append(neighbor_idx)

        if 2 <= chain_length <= 5:  # Short-chain (2 to 5 carbons)
            return True, "CoA moiety and short-chain fatty acyl thioester linkage present"

    return False, "Does not satisfy short-chain fatty acyl criteria, or invalid linkages"