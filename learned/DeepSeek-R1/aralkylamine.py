"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
An alkylamine in which the alkyl group is substituted by an aromatic group.
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine has an amine group connected to an alkyl chain which is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Find all amine groups (primary, secondary, tertiary; exclude amides)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    if not amine_matches:
        return False, "No amine group found"
    
    # Check each amine for connection to aromatic via alkyl chain (at least one carbon)
    for amine_match in amine_matches:
        n_idx = amine_match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        for neighbor in n_atom.GetNeighbors():
            # Only follow carbon chains
            if neighbor.GetAtomicNum() != 6:
                continue
                
            # Traverse the alkyl chain starting from the neighbor atom
            visited = set()
            stack = [(neighbor, 1)]  # (atom, depth from amine)
            
            while stack:
                current_atom, depth = stack.pop()
                if current_atom.GetIdx() in visited:
                    continue
                visited.add(current_atom.GetIdx())
                
                # Check if current atom is part of an aromatic system with sufficient depth
                if current_atom.GetIsAromatic() and depth >= 1:
                    return True, "Amine connected to aromatic group via alkyl chain"
                
                # Continue traversal through carbon atoms only
                if current_atom.GetAtomicNum() == 6:
                    for next_neighbor in current_atom.GetNeighbors():
                        # Avoid backtracking to previous atoms in the chain
                        if next_neighbor.GetIdx() != n_idx and next_neighbor.GetAtomicNum() == 6:
                            stack.append((next_neighbor, depth + 1))
    
    return False, "No aromatic group connected via alkyl chain to amine"