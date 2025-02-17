"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: decanoate ester
Definition: A fatty acid ester resulting from the condensation of decanoic acid (capric acid)
with an alcohol or phenol. Decanoic acid has 10 carbons (CH3(CH2)8COOH), so the acyl part 
in a decanoate ester should be exactly CH3(CH2)8CO–. 

Our approach:
1. Define two SMARTS patterns (to catch both orientations) for the ester fragment:
   Pattern A: "OC(=O)CCCCCCCCC" where the acyl part is on the right.
   Pattern B: "CCCCCCCCC(=O)O" where the acyl part is on the left.
2. For each match we check that the decanoate fragment is not part of a larger chain.
   In Pattern A the matched atom ordering is:
       index 0: the ester oxygen,
       index 1: the carbonyl carbon,
       index 2 to 10: the 9 alkyl carbons (with the first of these attached to the carbonyl)
   For a valid decanoate fragment, the terminal alkyl carbon (index 10) should be “terminal” 
   (i.e. have only a single carbon neighbor, coming from the chain) so that it is not part of a longer chain.
3. Similarly, for Pattern B the ordering is:
       index 0 to 8: the 9 alkyl carbons,
       index 9: the carbonyl carbon,
       index 10: the ester oxygen.
   In that case the terminal alkyl carbon is at index 0.
4. If any match passes the check we accept that the decanoate ester exists.
5. If none passes, then we return False.
   
This additional verification should help reduce false positives.
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester must contain an ester linkage in which the acyl (decanoate) portion 
    comes from decanoic acid (capric acid: CH3(CH2)8CO–, i.e. 10 carbons total).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a decanoate ester group with the proper chain length,
              False otherwise.
        str: Reason for the classification.
    """
    # Parse the input SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for the ester fragment.
    # Note: decanoate acyl group (CH3(CH2)8CO–) has a carbonyl carbon plus 9 alkyl carbons.
    # Pattern A: Ester oxygen on the left: O-C(=O)-CCCCCCCCC, where "CCCCCCCCC" is exactly 9 C's.
    patternA = Chem.MolFromSmarts("OC(=O)CCCCCCCCC")
    # Pattern B: Ester oxygen on the right: CCCCCCCCC(=O)O
    patternB = Chem.MolFromSmarts("CCCCCCCCC(=O)O")
    
    # Check matches for Pattern A.
    matches = mol.GetSubstructMatches(patternA)
    for match in matches:
        # match is a tuple of atom indices corresponding to the pattern.
        # In pattern A: index 0 -> ester O; index 1 -> carbonyl C; indices 2 to 10 -> 9 alkyl carbons.
        if len(match) != 11:
            continue  # not expected, but skip if length is off
        terminal_idx = match[10]  # the last alkyl carbon in the decanoate group
        terminal_atom = mol.GetAtomWithIdx(terminal_idx)
        # Check if terminal atom is really terminal (only bonded to the previous chain atom in the match).
        # Count neighbors that are carbon atoms.
        carbon_neighbors = [nbr for nbr in terminal_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # For a terminal CH3 normally it should be bonded only to one carbon (the previous in chain).
        if len(carbon_neighbors) == 1:
            return True, "Contains a decanoate ester functional group (matched via pattern A)"
    # Check matches for Pattern B.
    matches = mol.GetSubstructMatches(patternB)
    for match in matches:
        # In pattern B: indices 0-8 are the 9 alkyl carbons; index 9 is the carbonyl carbon; index 10 is the ester oxygen.
        if len(match) != 11:
            continue
        terminal_idx = match[0]  # here the terminal alkyl carbon is at the beginning of the chain
        terminal_atom = mol.GetAtomWithIdx(terminal_idx)
        # Check if the terminal alkyl carbon is only connected to one carbon (the next in the chain).
        carbon_neighbors = [nbr for nbr in terminal_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            return True, "Contains a decanoate ester functional group (matched via pattern B)"
    
    # If no valid decanoate ester fragment is found, return False.
    return False, "No valid decanoate ester group found in the molecule"

# Example usage (for testing purposes):
if __name__ == '__main__':
    # Test a set of example SMILES strings.
    test_smiles_list = [
        # Ethyl decanoate (should be true)
        ("O(C(CCCCCCCCC)=O)CC", "ethyl decanoate"),
        # Isobutyl decanoate (should be true)
        ("O(C(CCCCCCCCC)=O)CC(C)C", "isobutyl decanoate"),
        # A false positive example (one of the false positives given earlier) abbreviated:
        ("P(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCC)CO/C=C\\CCCCCCCCCCCCCC)(OCCN)(O)=O", "PE(P-16:0/15:1)"),
    ]
    for smi, name in test_smiles_list:
        result, reason = is_decanoate_ester(smi)
        print(f"Test Molecule: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")