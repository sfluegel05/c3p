"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:134741 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is formed by the condensation of decanoic acid with an alcohol/phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # SMARTS pattern for ester group: [OX2][CX3](=[OX1])
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if not ester_matches:
        return False, "No ester group found"
    
    # Check each ester group for decanoate (10-carbon chain on the acid side)
    for match in ester_matches:
        # Get the carbonyl carbon (index 1 in the SMARTS pattern)
        carbonyl_carbon = match[1]
        
        # Traverse the chain attached to the carbonyl carbon
        # Exclude the oxygen and the rest of the ester
        chain = []
        stack = [(carbonyl_carbon, 0)]  # (atom index, chain length)
        visited = set()
        
        while stack:
            atom_idx, length = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            
            # Stop at branching or non-carbon atoms
            if atom.GetDegree() > 2 or atom.GetAtomicNum() != 6:
                continue
            
            # Add neighbors (excluding the carbonyl oxygen)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx == match[0]:  # Skip the ester oxygen
                    continue
                if neighbor.GetAtomicNum() == 6:
                    stack.append((neighbor_idx, length + 1))
            
            chain.append(length)
        
        max_chain_length = max(chain) if chain else 0
        
        # Check if the chain length is 9 (since the carbonyl is the 10th carbon)
        if max_chain_length >= 9:
            return True, "Contains a decanoate ester group"
    
    return False, "No decanoate chain found in ester groups"