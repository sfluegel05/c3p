"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: decanoate ester
Definition: A fatty acid ester resulting from the formal condensation of the carboxy group 
of decanoic acid (capric acid, CH3(CH2)8COOH) with an alcohol or phenol. The decanoate acyl group 
must consist of exactly 10 carbons in length (CH3(CH2)8CO–) attached via an ester linkage.
Our improved approach:
1. Define two SMARTS patterns for the ester fragment:
   Pattern A: "OC(=O)CCCCCCCCC" (ester oxygen on left).
   Pattern B: "CCCCCCCCC(=O)O" (ester oxygen on right).
2. For each match (which returns 11 atoms) we check that the terminal alkyl carbon 
   (index 10 in pattern A, index 0 in pattern B) is terminal.
   • We ensure that all its carbon neighbors are part of the match.
   • A valid terminal CH3 should have exactly one carbon neighbor (the previous CH2 in the chain).
3. If a match passes these extra checks, we classify the molecule as containing a decanoate ester.
4. If no match passes, we return that no valid decanoate ester group exists.
"""

from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester must contain an ester linkage where the acyl (decanoate) portion 
    is exactly CH3(CH2)8CO–.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a decanoate ester group with the proper chain length,
              False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for the decanoate ester fragment.
    # A decanoate group (from decanoic acid) should have exactly 10 carbons:
    # CH3-(CH2)8-C(=O)-
    # Pattern A: ester oxygen is on the left: O-C(=O)-CCCCCCCCC
    patternA = Chem.MolFromSmarts("OC(=O)CCCCCCCCC")
    # Pattern B: ester oxygen is on the right: CCCCCCCCC(=O)O
    patternB = Chem.MolFromSmarts("CCCCCCCCC(=O)O")
    
    # Helper function to check if the terminal carbon (an alkyl CH3 group) is not part of a larger chain.
    def terminal_check(atom, match_indices):
        """
        Check that the given atom (terminal carbon in the candidate acyl chain) has exactly one carbon 
        neighbor and that any carbon neighbor it has is within the current match. This ensures the CH3 group 
        is not further substituted.
        """
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # Count only the neighbors that are in the match.
        in_match = [nbr for nbr in carbon_neighbors if nbr.GetIdx() in match_indices]
        # The terminal CH3 should have exactly one carbon neighbor (from the acyl chain) 
        # and no carbon neighbor outside the match.
        if len(in_match) == 1 and len(carbon_neighbors) == 1:
            return True
        return False

    # Check candidates for Pattern A.
    matches = mol.GetSubstructMatches(patternA)
    for match in matches:
        # For pattern A we expect 11 atoms: index 0: ester oxygen, index 1: carbonyl C, indices 2-10: 9 alkyl carbons.
        if len(match) != 11:
            continue  # unexpected match length
        terminal_idx = match[10]  # terminal alkyl carbon
        terminal_atom = mol.GetAtomWithIdx(terminal_idx)
        if terminal_check(terminal_atom, match):
            return True, "Contains a decanoate ester functional group (matched via pattern A)"
    
    # Check candidates for Pattern B.
    matches = mol.GetSubstructMatches(patternB)
    for match in matches:
        # For pattern B, indices 0-8 are the 9 alkyl carbons; index 9 is the carbonyl C; index 10 is ester oxygen.
        if len(match) != 11:
            continue
        terminal_idx = match[0]  # terminal alkyl carbon is at the beginning
        terminal_atom = mol.GetAtomWithIdx(terminal_idx)
        if terminal_check(terminal_atom, match):
            return True, "Contains a decanoate ester functional group (matched via pattern B)"
    
    # If no valid decanoate ester fragment is found, return False.
    return False, "No valid decanoate ester group found in the molecule"

# Example usage:
if __name__ == '__main__':
    # A few example SMILES (some expected to be true for decanoate esters and others as tests)
    test_smiles_list = [
        ("O(C(CCCCCCCCC)=O)CC", "ethyl decanoate"),
        ("O(C(CCCCCCCCC)=O)CC(C)C", "isobutyl decanoate"),
        ("CCCCCCCCCC(=O)OC", "methyl decanoate"),
        # A candidate that may be part of a longer chain (should be rejected)
        ("CCCCCCCCCCCC(=O)OC", "Undecyl ester (not decanoate)"),
    ]
    for smi, name in test_smiles_list:
        result, reason = is_decanoate_ester(smi)
        print(f"Test Molecule: {name}")
        print(f"SMILES: {smi}")
        print(f"Result: {result}")
        print(f"Reason: {reason}\n")