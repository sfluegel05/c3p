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
    
    # Check each ester group for decanoate (10 carbons in acid part including carbonyl)
    for match in ester_matches:
        carbonyl_carbon = match[1]
        visited = set([match[0]])  # Exclude ester oxygen from traversal
        stack = [carbonyl_carbon]
        carbon_count = 0
        
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue  # Only count carbons
            
            carbon_count += 1
            
            # Add all adjacent atoms except ester oxygen
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited and neighbor.GetIdx() != match[0]:
                    stack.append(neighbor.GetIdx())
        
        # Need at least 10 carbons in acid chain (carbonyl + 9 more)
        if carbon_count >= 10:
            return True, "Contains decanoate chain in ester group"
    
    return False, "No ester group with 10-carbon acid chain found"