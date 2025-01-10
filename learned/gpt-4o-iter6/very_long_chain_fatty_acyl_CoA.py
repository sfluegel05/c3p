"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Classifies a molecule as a very long-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Convert the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Substructure pattern for Coenzyme A (CoA) part, using a comprehensive pattern
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)O[C@H]1[C@H](O)[C@H](O)C(CO1)OP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"
    
    # Substructure pattern for the thioester linkage: carbonyl group connected to a sulfur atom
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)S")
    carbonyl_match = mol.GetSubstructMatches(carbonyl_pattern)
    if not carbonyl_match:
        return False, "No carbonyl linked to CoA structure found"

    # Track longest consecutive chain of carbon atoms from each carbonyl carbon
    max_chain_length = 0
    visited = set()

    # Define depth-first search function to track acyl chain
    def find_chain_length(atom_idx, length):
        nonlocal max_chain_length
        if length > max_chain_length:
            max_chain_length = length
        visited.add(atom_idx)
        
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:  # Check for carbon atoms
                find_chain_length(neighbor.GetIdx(), length + 1)

    # Loop over all carbonyl matches to explore different potential chains
    for match in carbonyl_match:
        carbonyl_carbon_idx = match[0]
        visited.clear()  # Clear visited sets for each new carbonyl group
        for neighbor in mol.GetAtomWithIdx(carbonyl_carbon_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:  # Start DFS from carbon atom
                find_chain_length(neighbor.GetIdx(), 1)

    # Check if the longest carbon chain exceeds C22
    if max_chain_length <= 22:
        return False, f"Longest carbon chain is {max_chain_length}, must be greater than 22"

    return True, "Contains CoA structure with a fatty acyl chain longer than C22"