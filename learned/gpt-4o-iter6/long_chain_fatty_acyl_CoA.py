"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA includes a long-chain fatty acid esterified with the thiol group of coenzyme A.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Redefine Coenzyme A motif with a more comprehensive pattern
    coa_pattern = Chem.MolFromSmarts("C(=O)NCCC(=O)NCCSC(=O)CCNC(=O)[C@H](O)[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc23)[C@H](O)[C@H]1O)C(C)(C)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Define thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"
    
    # Check the fatty acid part for carbon chain length
    valid_long_chain = False
    for match in thioester_matches:
        # Assume the carbon following the thioester is the beginning of the chain
        start_atom = match[2]  # Index of the first carbon in the chain
        
        # Count the number of carbons in the fatty acid chain
        carbon_chain = set()
        queue = [mol.GetAtomWithIdx(start_atom)]
        visited = set([start_atom])

        while queue:
            current_atom = queue.pop()
            # Consider all continuous carbon atoms
            if current_atom.GetSymbol() == 'C':
                carbon_chain.add(current_atom.GetIdx())

            for neighbor in current_atom.GetNeighbors():
                if neighbor.GetIdx() not in visited and neighbor.GetSymbol() == 'C':
                    visited.add(neighbor.GetIdx())
                    queue.append(neighbor)

        carbon_chain_length = len(carbon_chain)
        # Check if chain length is within the specified range
        if 13 <= carbon_chain_length <= 22:
            valid_long_chain = True
            break
            
    if not valid_long_chain:
        return False, f"Carbon chain length not in valid range (C13-C22), got {carbon_chain_length}"
    
    return True, f"Valid long-chain fatty acyl-CoA with {carbon_chain_length} carbon atoms"