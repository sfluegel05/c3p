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
    
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage required for acyl-CoA not found"
    
    # Estimate the fatty acyl chain length (13 to 22)
    # Assuming the attachment is through the thioester carbon
    thio_linkages = mol.GetSubstructMatches(thioester_pattern)
    for linkage in thio_linkages:
        carbon_atom = mol.GetAtomWithIdx(linkage[0])
        carbon_chain_length = 0
        visited = set()
        
        # We will perform a breadth-first search (BFS) to count the carbon chain
        atoms_to_visit = [carbon_atom]
        while atoms_to_visit:
            current_atom = atoms_to_visit.pop(0)
            if current_atom.GetIdx() not in visited and current_atom.GetAtomicNum() == 6:  # carbon check
                visited.add(current_atom.GetIdx())
                carbon_chain_length += 1
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited:
                        atoms_to_visit.append(neighbor)
        
        # Check if the carbon_chain_length is within the specified range (13 to 22)
        if 13 <= carbon_chain_length <= 22:
            return True, "Structure matches long-chain fatty acyl-CoA"
    
    return False, "Fatty acyl chain not in the correct length range (13-22 carbons)"