"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: decanoate ester

A decanoate ester is defined as a fatty acid ester resulting from the condensation 
of decanoic acid (capric acid: CH3(CH2)8COOH) with an alcohol or phenol.
In the acyl group the decanoyl fragment must be exactly:
    CH3-(CH2)8-C(=O)-
We search for two SMARTS patterns:
    Pattern A: ester oxygen on the left -> "OC(=O)CCCCCCCCC"
    Pattern B: ester oxygen on the right -> "CCCCCCCCC(=O)O"
Then, for each candidate we verify:
    - The candidate match returns 11 atoms.
    - In mode A:
         index0: O (ester oxygen) [may be linked externally]
         index1: carbonyl C (should connect to index0 and index2)
         index2-9: chain methylenes (each connects to previous and next; must not show extra carbon branching)
         index10: terminal CH3 (must only connect to index9 within the chain)
    - In mode B the order is reversed.
If any candidate passes these strict tests, we return a positive classification.
"""

from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    The decanoate moiety (decanoyl fragment) must appear exactly as:
        CH3-(CH2)8-C(=O)-  (i.e. 10 carbons total, with a terminal methyl and no branching)
    and be connected via an ester linkage.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a valid decanoate ester functional group, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for the decanoate acyl chain.
    # Pattern A: the ester O appears before the carbonyl:
    #    [O]-[C](=O)-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]
    patternA = Chem.MolFromSmarts("OC(=O)CCCCCCCCC")
    # Pattern B: the ester O appears after the acyl chain:
    #    [C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C](=O)-[O]
    patternB = Chem.MolFromSmarts("CCCCCCCCC(=O)O")
    
    # Helper: check connectivity and bonding in a candidate decanoate match.
    def check_candidate(match, mode):
        # Expect match length == 11 atoms in both patterns.
        if len(match) != 11:
            return False
        
        # In order to simplify checking, we retrieve the atoms in order.
        # For mode "A": indices:
        #   0: ester O, 1: carbonyl C, 2-10: decanoate chain (10 carbons, with index10 being terminal CH3)
        # For mode "B": indices:
        #   0: terminal CH3, 1-8: chain methylenes, 9: carbonyl C, 10: ester O
        # We will verify element types and connectivity along the chain.
        if mode == "A":
            # Check element types:
            if mol.GetAtomWithIdx(match[0]).GetSymbol() != "O":
                return False
            if mol.GetAtomWithIdx(match[1]).GetSymbol() != "C":
                return False
            for i in range(2, 11):
                if mol.GetAtomWithIdx(match[i]).GetSymbol() != "C":
                    return False
            # Connectivity along the chain:
            # The ester oxygen (index0) must be bonded to the carbonyl carbon (index1).
            if not mol.GetBondBetweenAtoms(match[0], match[1]):
                return False
            # The carbonyl carbon (index1) should bond to index0 and index2.
            if not (mol.GetBondBetweenAtoms(match[1], match[2]) and mol.GetBondBetweenAtoms(match[1], match[0])):
                return False
            # Check the linear chain from index2 to index10.
            # For indices 2 to 9 (middle chain atoms): each should have bonds to the previous and next atom (in the match)
            for i in range(2, 10):
                if not (mol.GetBondBetweenAtoms(match[i], match[i-1]) and mol.GetBondBetweenAtoms(match[i], match[i+1])):
                    return False
                # Also ensure that these atoms do not have extra bonds to carbons that are also in the match
                for nbr in mol.GetAtomWithIdx(match[i]).GetNeighbors():
                    if nbr.GetSymbol() == "C" and nbr.GetIdx() in match:
                        if nbr.GetIdx() not in (match[i-1], match[i+1]):
                            return False
            # Terminal CH3 (index 10): should bond only to index9 among atoms in the match.
            if not mol.GetBondBetweenAtoms(match[10], match[9]):
                return False
            # Check that terminal CH3 is truly terminal (only one C neighbor in the molecule that is in the match).
            terminal_atom = mol.GetAtomWithIdx(match[10])
            c_neighbors_in_match = [nbr.GetIdx() for nbr in terminal_atom.GetNeighbors() if nbr.GetSymbol() == "C" and nbr.GetIdx() in match]
            if len(c_neighbors_in_match) != 1:
                return False
            # Optionally, check that the terminal CH3 has exactly 3 attached hydrogens.
            if terminal_atom.GetTotalNumHs() != 3:
                return False
            return True
        
        elif mode == "B":
            # For mode B, expected order:
            #  index0: terminal CH3, indices1-8: chain methylenes, index9: carbonyl C, index10: ester O.
            if mol.GetAtomWithIdx(match[10]).GetSymbol() != "O":
                return False
            if mol.GetAtomWithIdx(match[9]).GetSymbol() != "C":
                return False
            for i in range(0, 9):
                if mol.GetAtomWithIdx(match[i]).GetSymbol() != "C":
                    return False
            # Connectivity:
            # Terminal CH3 (index0) should bond only to index1.
            if not mol.GetBondBetweenAtoms(match[0], match[1]):
                return False
            terminal_atom = mol.GetAtomWithIdx(match[0])
            c_neighbors_in_match = [nbr.GetIdx() for nbr in terminal_atom.GetNeighbors() if nbr.GetSymbol() == "C" and nbr.GetIdx() in match]
            if len(c_neighbors_in_match) != 1:
                return False
            if terminal_atom.GetTotalNumHs() != 3:
                return False
            # For indices 1 to 8 (chain methylenes): each should bond to the previous and next.
            for i in range(1, 9):
                if not (mol.GetBondBetweenAtoms(match[i], match[i-1]) and mol.GetBondBetweenAtoms(match[i], match[i+1])):
                    return False
                # Ensure no extra in-chain carbon branches.
                for nbr in mol.GetAtomWithIdx(match[i]).GetNeighbors():
                    if nbr.GetSymbol() == "C" and nbr.GetIdx() in match:
                        if nbr.GetIdx() not in (match[i-1], match[i+1]):
                            return False
            # Carbonyl carbon at index9 must bond to index8 and index10.
            if not (mol.GetBondBetweenAtoms(match[9], match[8]) and mol.GetBondBetweenAtoms(match[9], match[10])):
                return False
            return True
        
        return False

    # Check candidate matches for Pattern A.
    matchesA = mol.GetSubstructMatches(patternA)
    for match in matchesA:
        if check_candidate(match, "A"):
            return True, "Contains a decanoate ester functional group (matched via pattern A)"
    
    # Check candidate matches for Pattern B.
    matchesB = mol.GetSubstructMatches(patternB)
    for match in matchesB:
        if check_candidate(match, "B"):
            return True, "Contains a decanoate ester functional group (matched via pattern B)"
    
    return False, "No valid decanoate ester group found in the molecule"

# Example usage:
if __name__ == '__main__':
    # Test a few examples from the provided list.
    examples = [
        ("O(C(CCCCCCCCC)=O)CC", "ethyl decanoate"),
        ("O(C(CCCCCCCCC)=O)CC(C)C", "isobutyl decanoate"),
        ("O=C(OC)CCCCCCCCC", "methyl decanoate"),
        ("CCCCCCCCC(=O)OCC", "alternative ethyl decanoate"),
        # More complex examples:
        ("CCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCC", "1-palmitoyl-2-decanoyl-sn-glycero-3-phosphocholine"),
        ("CCCCCCCCC(=O)OC[C@@H](O)CO", "1-decanoyl-sn-glycerol"),
    ]
    for smi, name in examples:
        res, reason = is_decanoate_ester(smi)
        print(f"Test Molecule: {name}")
        print(f"SMILES: {smi}")
        print(f"Result: {res}, Reason: {reason}\n")