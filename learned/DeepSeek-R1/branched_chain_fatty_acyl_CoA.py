"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:22882 branched-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA has a CoA moiety attached via a thioester bond to a branched fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define corrected thioester pattern (S-C(=O)-C without valence restriction)
    thioester_pattern = Chem.MolFromSmarts("[S]-C(=O)-C")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found (not an acyl-CoA)"
    
    # Check each thioester match for branched R group
    for match in thioester_matches:
        if len(match) != 3:
            continue
        
        # Get indices of S, C=O, and first R-group carbon
        s_idx, co_idx, r_start_idx = match
        
        # Traverse R-group (fatty acid chain) starting from r_start_idx
        visited = {s_idx, co_idx}
        stack = [r_start_idx]
        r_group_atoms = set()
        
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            r_group_atoms.add(atom_idx)
            
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited:
                    stack.append(neighbor.GetIdx())
        
        # Check for branching in R-group (any carbon with degree >=3)
        for atom_idx in r_group_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6 and atom.GetDegree() >= 3:
                return True, "Branched chain detected in fatty acid group"
    
    return False, "No branching detected in fatty acid chain"