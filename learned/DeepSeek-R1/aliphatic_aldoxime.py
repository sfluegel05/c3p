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
    An aliphatic aldoxime is derived from an aliphatic aldehyde, containing R-CH=N-OH where R is an aliphatic, acyclic group.

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

    # Get ring information for later checks
    ring_info = mol.GetRingInfo()
    ring_info.Initialize()  # Ensure ring info is calculated

    # Check each aldoxime group for aliphatic characteristics
    for match in aldoxime_matches:
        cx3_idx = match[0]  # Aldoxime carbon (C in C=N-OH)
        cx3_atom = mol.GetAtomWithIdx(cx3_idx)
        
        # Find non-hydrogen, non-nitrogen neighbors (must be exactly 1 for aldehyde-derived)
        r_candidates = [n for n in cx3_atom.GetNeighbors()
                       if n.GetAtomicNum() == 6 and n.GetIdx() != match[1]]
        
        # Check for exactly one carbon neighbor (R group)
        if len(r_candidates) != 1:
            continue
            
        r_carbon = r_candidates[0]
        
        # Check R group is acyclic and aliphatic
        visited = set()
        stack = [r_carbon]
        is_valid = True
        
        while stack and is_valid:
            current = stack.pop()
            
            # Check for aromatic atoms or rings in R group
            if current.GetIsAromatic() or ring_info.IsAtomInRing(current.GetIdx()):
                is_valid = False
                break
            
            # Traverse all connected atoms except aldoxime carbon
            for neighbor in current.GetNeighbors():
                if neighbor.GetIdx() == cx3_idx:
                    continue
                if neighbor.GetIdx() not in visited:
                    visited.add(neighbor.GetIdx())
                    stack.append(neighbor)
        
        if is_valid:
            # Verify aldoxime carbon has only one non-H neighbor (the R group)
            non_h_neighbors = [n for n in cx3_atom.GetNeighbors()
                              if n.GetAtomicNum() != 1 and n.GetIdx() != match[1]]
            if len(non_h_neighbors) == 1:
                return True, "Contains aldoxime group with acyclic aliphatic R chain"

    return False, "No valid aliphatic aldoxime groups found"