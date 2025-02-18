"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: CHEBI:83406 aliphatic aldoxime
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

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

    # Check each aldoxime group for aliphatic characteristics
    for match in aldoxime_matches:
        cx3_idx = match[0]  # Aldoxime carbon (C in C=N-OH)
        cx3_atom = mol.GetAtomWithIdx(cx3_idx)
        
        # Find non-hydrogen, non-nitrogen neighbors (must be exactly 1 for aldehyde-derived)
        r_candidates = [n for n in cx3_atom.GetNeighbors()
                       if n.GetAtomicNum() != 1 and n.GetIdx() != match[1]]  # Exclude H and N
        
        # Check for exactly one R group neighbor (should be carbon)
        if len(r_candidates) != 1:
            continue
            
        r_group = r_candidates[0]
        if r_group.GetAtomicNum() != 6:
            continue  # R group must be carbon-based
        
        # Check R group is aliphatic (no aromatic bonds) and acyclic (no rings)
        # Traverse all atoms in R group starting from connection point
        visited = set()
        stack = [r_group]
        is_aliphatic = True
        has_ring = False
        
        while stack and is_aliphatic and not has_ring:
            current = stack.pop()
            if current.GetIdx() in visited:
                continue
            visited.add(current.GetIdx())
            
            # Check for aromatic atoms or bonds
            if current.GetIsAromatic():
                is_aliphatic = False
                break
            
            # Check bonds for aromaticity
            for bond in current.GetBonds():
                if bond.GetBondType() == Chem.BondType.AROMATIC:
                    is_aliphatic = False
                    break
            if not is_aliphatic:
                break
            
            # Check if current atom is in a ring
            if mol.GetRingInfo().IsAtomInRing(current.GetIdx()):
                has_ring = True
                break
            
            # Add neighboring atoms except aldoxime carbon
            for neighbor in current.GetNeighbors():
                if neighbor.GetIdx() != cx3_idx and neighbor.GetIdx() not in visited:
                    stack.append(neighbor)
        
        if has_ring or not is_aliphatic:
            continue
        
        # Check that aldoxime carbon has no other substituents (only R and N=O)
        # All other bonds should be to hydrogen
        h_count = 0
        for bond in cx3_atom.GetBonds():
            other_atom = bond.GetOtherAtom(cx3_atom)
            if other_atom.GetAtomicNum() == 1:
                h_count += 1
            elif other_atom.GetIdx() not in [r_group.GetIdx(), match[1]]:
                # Extra substituent found (not H, R, or N)
                break
        else:
            # All checks passed
            return True, "Contains aldoxime group with acyclic aliphatic R chain"

    return False, "No valid aliphatic aldoxime groups found"