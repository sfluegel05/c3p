"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.

    A short-chain fatty acyl-CoA is a fatty acyl-CoA formed by the condensation
    of the thiol group of coenzyme A with a short-chain fatty acid.

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
    
    # A more flexible search for the CoA moiety based on common substructure without assuming full structure
    # We focus only on crucial portions like the pantetheine and the phosphoadenosine diphosphate
    coa_pantetheine_pattern = Chem.MolFromSmarts("NCC(=O)NCCSC")
    coa_diphosphate_pattern = Chem.MolFromSmarts("O=P(O)(O)OP(O)(O)=O")

    if not (mol.HasSubstructMatch(coa_pantetheine_pattern) and mol.HasSubstructMatch(coa_diphosphate_pattern)):
        return False, "No definitive CoA moiety found"
    
    # Detect thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    
    if not thioester_matches:
        return False, "No thioester linkage found"
    
    # Ensure short-chain fatty acid component (2-6 carbons)
    for match in thioester_matches:
        sulfur = mol.GetAtomWithIdx(match[1])
        connected_carbons = [atom for atom in sulfur.GetNeighbors() if atom.GetSymbol() == 'C']
        if not connected_carbons:
            continue
        
        # Start the carbon chain count from the carbon bonded to sulfur
        chain_carbon = connected_carbons[0]
        
        # Use a BFS or DFS to count carbons in the chain
        visited = set()
        carbon_count = 0
        stack = [chain_carbon]
        
        while stack:
            atom = stack.pop()
            if atom.GetIdx() not in visited:
                visited.add(atom.GetIdx())
                if atom.GetSymbol() == 'C':
                    carbon_count += 1
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in visited:
                            stack.append(neighbor)
        
        if 2 <= carbon_count <= 6:
            return True, "Contains CoA moiety and a compatible short-chain fatty acyl group"
    
    return False, "No proper short-chain fatty acyl component found"