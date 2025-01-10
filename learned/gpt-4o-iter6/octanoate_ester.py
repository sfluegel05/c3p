"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester features the presence of the octanoyl group, which is a C(=O)O group followed by an 8-carbon chain.

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

    # SMARTS pattern for ester linkages (generic)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(ester_pattern)
    
    if not matches:
        return False, "No ester linkage found"
    
    # Analyze matches for octanoyl characteristics
    for match in matches:
        # Get the atom index after the ester oxygen to count carbon chain length
        oxygen_idx = match[2]
        carbon_chain_length = 0
        visited = set()
        to_visit = [oxygen_idx]
        
        while to_visit:
            curr_idx = to_visit.pop()
            if curr_idx in visited:
                continue
            visited.add(curr_idx)

            # Get neighbors
            curr_atom = mol.GetAtomWithIdx(curr_idx)
            if curr_atom.GetAtomicNum() == 6:  # It's a carbon
                carbon_chain_length += 1
            for neighbor in curr_atom.GetNeighbors():
                neigh_idx = neighbor.GetIdx()
                if neighbor.GetAtomicNum() == 6 and neigh_idx not in visited:
                    to_visit.append(neigh_idx)
        
        if carbon_chain_length == 8:
            return True, "Contains valid octanoyl ester linkage with 8-carbon chain"
    
    return False, "No valid octanoyl ester linkage with an 8-carbon chain found"