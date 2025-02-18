"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: Octanoate ester 
Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid)
i.e. the acyl portion is CH3(CH2)6C(=O)O–.
"""

from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    
    An octanoate ester is defined as any fatty acid ester where the acyl group 
    originates from octanoic acid (CH3(CH2)6C(=O)O–). 
    We first search for the SMARTS pattern representing this group. Then, for each match,
    we enforce that the matched acyl chain has the expected chain length and is terminal
    (i.e. the terminal methyl does not have bonds to extra carbons).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a proper octanoate ester substructure is present, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the basic SMARTS for an octanoate ester acyl group.
    # Octanoate (acyl) group: CH3(CH2)6C(=O)O ; note that the SMARTS "CCCCCCCC(=O)O" 
    # matches 8 contiguous C atoms, with the 8th functioning as the carbonyl (C(=O)) attached to O.
    octanoate_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O")
    if octanoate_pattern is None:
        return False, "Internal error: invalid SMARTS pattern"
    
    # Get all substructure matches.
    matches = mol.GetSubstructMatches(octanoate_pattern)
    if not matches:
        return False, "No octanoate ester group found"
    
    # For each match, do extra validation of the acyl chain.
    # We expect the match to map the acyl chain where:
    #  - The first atom (index 0 in the match) should be a terminal CH3 (only one heavy neighbor).
    #  - The carbonyl carbon (expected to be at index 7 in the SMARTS, 
    #    since "CCCCCCCC(=O)O" gives eight C atoms followed by the ester oxygen)
    #    should be bonded only to the preceding carbon (index 6) and to the ester O.
    for match in matches:
        # In our SMARTS, match[0] is the terminal methyl of the acyl chain.
        terminal_idx = match[0]
        terminal_atom = mol.GetAtomWithIdx(terminal_idx)
        # Check that the terminal CH3 is not connected to any extra heavy atom outside the match.
        heavy_neighbors = [nbr.GetIdx() for nbr in terminal_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # We expect only one neighbor and it should be part of the match.
        if len(heavy_neighbors) != 1 or heavy_neighbors[0] not in match:
            continue  # not a proper terminal group
        
        # Next, check the carbonyl carbon.
        # In the pattern, the carbonyl carbon is the 8th atom in the chain.
        # Note: Depending on implicit handling, the match tuple length should be 9 
        # (8 carbons and the ester oxygen). We check if we have at least 8 mapped carbons.
        if len(match) < 8:
            continue  # unexpected mapping
        carbonyl_idx = match[7]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Get neighbors’ indices (only consider heavy atoms).
        carbonyl_neighbors = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # We expect the carbonyl carbon to be bonded to exactly 2 atoms:
        # one from the acyl chain (should be match[6]) and one ester oxygen (should be part of the match).
        if len(carbonyl_neighbors) != 2 or (match[6] not in carbonyl_neighbors):
            continue
        
        # If we passed both tests, we consider the match a valid octanoate ester.
        return True, "Contains an octanoate ester group (acyl derived from octanoic acid) with proper chain length and terminal methyl."
    
    # If none of the matches qualifies for extra constraints, return false.
    return False, "Octanoate ester-like substructure found, but chain length/terminal conditions not met."