"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is characterized by the ester group (-C(=O)O-) where the acid part 
    is specifically octanoic acid (8 carbon chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ester group pattern (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester linkage found"
    
    # For each ester match, check for 8-carbon chain attached to the carbonyl
    for match in ester_matches:
        # Identify the atom indices
        carbonyl_c, o_ester = match[0], match[2]
        
        # Traverse the chain from the ester oxygen (should lead back to carbonyl carbon)
        # Check the length and linearity of the carbon chain attached to the carbonyl
        chain_length = 0
        visited = set()
        to_visit = [o_ester]
        
        while to_visit:
            atom_idx = to_visit.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                chain_length += 1
                # Add connected atoms excluding those already visited or back to carbonyl C
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor_idx != carbonyl_c:
                        to_visit.append(neighbor_idx)
        
        # For octanoic acid, we should find exactly 8 carbons in a linear chain
        if chain_length == 8:
            return True, "Contains octanoate ester group"
    
    return False, "Ester group found, but not octanoate (8-carbon chain)"