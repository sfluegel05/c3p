"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: CHEBI:83406 aliphatic aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is derived from an aliphatic aldehyde, containing R-CH=N-OH where R is an aliphatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define aldoxime pattern: Carbon double bonded to nitrogen with hydroxyl group
    aldoxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2]-[OH]")
    if not aldoxime_pattern:
        return False, "Invalid pattern definition"

    aldoxime_matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not aldoxime_matches:
        return False, "No aldoxime group (C=N-OH) found"

    # Check each aldoxime group for aliphatic characteristics
    for match in aldoxime_matches:
        cx3_idx = match[0]  # Aldoxime carbon (C in C=N-OH)
        cx3_atom = mol.GetAtomWithIdx(cx3_idx)
        
        # Find adjacent carbon (R group in R-CH=N-OH)
        carbon_neighbors = [n for n in cx3_atom.GetNeighbors() 
                           if n.GetAtomicNum() == 6 and n.GetIdx() != match[1]]
        
        # Must have exactly one carbon neighbor (aldehyde-derived structure)
        if len(carbon_neighbors) != 1:
            continue
            
        r_carbon = carbon_neighbors[0]
        
        # Check if R group is aliphatic (non-aromatic and no aromatic bonds)
        visited = set()
        stack = [r_carbon]
        is_aliphatic = True
        
        while stack and is_aliphatic:
            current = stack.pop()
            if current.GetIsAromatic():
                is_aliphatic = False
                break
            # Traverse all connected atoms except aldoxime carbon
            for neighbor in current.GetNeighbors():
                if neighbor.GetIdx() == cx3_idx:
                    continue
                if neighbor.GetIdx() not in visited:
                    visited.add(neighbor.GetIdx())
                    stack.append(neighbor)
        
        if is_aliphatic:
            return True, "Contains aldoxime group with aliphatic R chain"

    return False, "No aliphatic aldoxime groups found"