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
    if not coa_pattern:
        return False, "Error creating CoA SMARTS pattern"
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A (CoA) moiety not found"
    
    # Define SMARTS pattern for the thioester linkage indicative of fatty acyl-CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not thioester_pattern:
        return False, "Error creating thioester SMARTS pattern"
    
    thio_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thio_matches:
        return False, "Thioester linkage required for acyl-CoA not found"
    
    # Estimate the longest unbranched fatty acyl chain length; assume attachment through thioester carbon
    for match in thio_matches:
        c_atom_idx = match[0]  # Carbon attached to thioester
        c_atom = mol.GetAtomWithIdx(c_atom_idx)
        
        acyl_chain_atoms = set()
        visited = set([c_atom_idx])
        stack = [c_atom]

        # Enhanced depth-first search (DFS)
        while stack:
            current_atom = stack.pop()
            if current_atom.GetIdx() not in visited:
                visited.add(current_atom.GetIdx())
                if current_atom.GetAtomicNum() == 6:  # Carbon filter
                    acyl_chain_atoms.add(current_atom.GetIdx())
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                            stack.append(neighbor)
        
        acyl_chain_length = len(acyl_chain_atoms)
        
        # Check if the length of the linear chain is within the designated range (13 to 22 carbons)
        if 13 <= acyl_chain_length <= 22:
            return True, "Structure matches long-chain fatty acyl-CoA"
    
    return False, "Fatty acyl chain not in the correct length range (13-22 carbons)"