"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the formal condensation of the thiol group 
    of coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined Coenzyme A pattern: Check presence of key coenzyme A structure
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)O[C@H]1O[C@H](CO[P]([O-])([O-]))[C@H](O)[C@H]1OP(=O)([O-])[O-]")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A moiety found"
    
    # Thioester linkage pattern for recognizing fatty acyl attachment
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCN")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Detect carbon chain length, starting from thioester carbon (first atom in match)
    # Traverse the molecule to determine the length of the longest contiguous carbon chain
    for match in thioester_matches:
        start_atom_idx = match[0]  # C(=O) carbon atom index

        # Use a breadth-first search (BFS) to determine carbon chain length
        visited = set()
        queue = [(start_atom_idx, 0)]  # (atom index, chain length)
        max_chain_length = 0
        
        while queue:
            current_idx, current_length = queue.pop(0)
            if current_idx in visited:
                continue
            
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            
            if atom.GetAtomicNum() == 6:  # Check for carbon
                max_chain_length = max(max_chain_length, current_length)

            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    queue.append((neighbor_idx, current_length + (1 if neighbor.GetAtomicNum() == 6 else 0)))

        # Check if the detected longest carbon chain is in the range of 13 to 22
        if 13 <= max_chain_length <= 22:
            return True, f"Valid long-chain fatty acyl-CoA with {max_chain_length} carbon atoms in chain"

    return False, "No suitable long-chain fatty acid chain found"