"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    
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
    
    # Define SMARTS pattern for CoA moiety
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A (CoA) moiety not found"
    
    # Define SMARTS pattern for the thioester linkage indicative of fatty acyl-CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thio_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thio_matches:
        return False, "Thioester linkage required for acyl-CoA not found"
    
    # Estimate the longest unbranched fatty acyl chain length; assume attachment through thioester carbon
    max_chain_length = 0
    for match in thio_matches:
        c_atom_idx = match[0]  # Carbon attached to thioester
        c_atom = mol.GetAtomWithIdx(c_atom_idx)

        # Traverse to find longest carbon chain
        chain_length = 0
        
        visited = set()
        stack = [(c_atom_idx, 1)]  # Store atom index and current chain length

        while stack:
            current_idx, current_length = stack.pop()
            
            # We only consider carbon atoms for chain length
            current_atom = mol.GetAtomWithIdx(current_idx)
            if current_atom.GetAtomicNum() == 6:  # If Carbon
                chain_length = max(chain_length, current_length)
            
            # Visit neighbors extending the carbon chain
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                    visited.add(neighbor_idx)
                    stack.append((neighbor_idx, current_length + 1))

        # Keep track of the maximum length found
        max_chain_length = max(max_chain_length, chain_length)
        
    # Check if the longest acyl chain is within 13 to 22 carbons
    if 13 <= max_chain_length <= 22:
        return True, "Structure matches long-chain fatty acyl-CoA with chain length in range"
    
    return False, "Fatty acyl chain not in the correct length range (13-22 carbons)"