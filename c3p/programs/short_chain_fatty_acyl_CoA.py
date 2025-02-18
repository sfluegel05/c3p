"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.

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

    # Look for Coenzyme A specific fragments
    coA_pattern = Chem.MolFromSmarts("OP(=O)(O)OC[C@H]1O[C@H]([C@@H](O)[C@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Coenzyme A structure not found"

    # Look for thioester linkage: S-C(=O) for fatty acyl attachment
    thioester_pattern = Chem.MolFromSmarts("S-[CX3](=O)-")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Check if the chain length from the carbon in S-C(=O)-C linkage is 2 to 6 carbons
    for match in thioester_matches:
        carbonyl_carbon_idx = match[1]
        next_carbon_idx = match[2]

        # Use BFS to traverse the carbon chain starting from next_carbon
        queue = [next_carbon_idx]
        visited = {next_carbon_idx}
        chain_length = 1

        while queue:
            current_idx = queue.pop(0)
            current_atom = mol.GetAtomWithIdx(current_idx)

            for neighbor in current_atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                n_atom = mol.GetAtomWithIdx(n_idx)

                # Visit only carbon atoms that are part of the chain and not revisiting
                if n_idx not in visited and n_atom.GetAtomicNum() == 6:
                    visited.add(n_idx)
                    queue.append(n_idx)
                    chain_length += 1
        
        # Validate if chain length is 2 to 6
        if 2 <= chain_length <= 6:
            return True, "Contains Coenzyme A linked with a short-chain fatty acid"

    return False, "No short-chain fatty acid attached to CoA"