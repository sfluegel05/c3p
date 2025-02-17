"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: decanoate ester

Definition:
 A decanoate ester is defined as a fatty acid ester resulting from the formal condensation 
 of decanoic acid (capric acid: CH3(CH2)8COOH) with an alcohol or phenol.
 In the acyl group (decanoate) the fragment must be the decanoyl fragment:
  CH3-(CH2)8-C(=O)-    (exactly 10 carbons: a terminal CH3 followed by 8 CH2 groups and a carbonyl C).

Our improved approach:
 1. Define two SMARTS patterns:
    Pattern A (ester oxygen on the left): "OC(=O)CCCCCCCCC"
      Expected match order: [O, C, C, C, C, C, C, C, C, C, C]
         where atom0 = ester O,
               atom1 = carbonyl C,
               atoms2-9 = eight chain carbons (CH2),
               atom10 = terminal CH3.
    Pattern B (ester oxygen on the right): "CCCCCCCCC(=O)O"
      Expected match order: [C, C, C, C, C, C, C, C, C, C, O]
         where atom0 = terminal CH3,
               atoms1-8 = eight chain carbons,
               atom9 = carbonyl C,
               atom10 = ester O.
 2. For each match (if its length is exactly 11) we check:
      - That the atoms in the positions are of the expected element type.
      - That the candidate “chain” is internally connected (i.e. there is a bond between consecutive atoms in the match).
      - That the terminal CH3 atom is truly terminal (only one bonded carbon) and carries exactly 3 hydrogens.
 3. If any candidate passes these tests, we classify the molecule as containing a decanoate ester.
"""

from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    The decanoate moiety must appear as the acyl group: CH3-(CH2)8-C(=O)– and be attached via an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a valid decanoate ester functional group, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for the decanoate group.
    # Pattern A: ester oxygen on the left.
    patternA = Chem.MolFromSmarts("OC(=O)CCCCCCCCC")
    # Pattern B: ester oxygen on the right.
    patternB = Chem.MolFromSmarts("CCCCCCCCC(=O)O")
    
    # Helper function: check that a terminal carbon candidate is indeed a terminal CH3.
    def terminal_atom_check(atom, match_indices):
        # Must be carbon.
        if atom.GetSymbol() != "C":
            return False
        # Should have exactly one neighbor that is carbon and that neighbor must be in the match.
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            return False
        # Total hydrogens (implicit+explicit) should be exactly 3.
        if atom.GetTotalNumHs() != 3:
            return False
        # Also, the single carbon neighbor must be in the candidate substructure.
        if carbon_neighbors[0].GetIdx() not in match_indices:
            return False
        return True

    # Helper function: verify linear connectivity of the candidate match.
    # expected_order is a list of tuples: (expected element, is_terminal, expected_neighbor_counts)
    # For a proper decanoate chain the connectivity (within the candidate atoms) should be:
    # For Pattern A:
    #  index 0 (ester O): connected only to index 1.
    #  index 1 (carbonyl C): connected to index 0 and index 2.
    #  index 2..9 (chain carbons): each connected to the preceding and following (2 neighbors).
    #  index 10 (terminal CH3): connected only to index 9.
    #
    # For Pattern B the order is reversed:
    #  index 0 (terminal CH3): connected only to index 1.
    #  index 1..8 (chain carbons): connected to preceding and following.
    #  index 9 (carbonyl C): connected to index 8 and index 10.
    #  index 10 (ester O): connected only to index 9.
    def check_candidate(match, pattern_type):
        if len(match) != 11:
            return False
        # Create a set for quick look-up
        match_set = set(match)
        # For each consecutive pair in the expected chain, verify there is a bond between them.
        if pattern_type == "A":
            # Expected order for Pattern A:
            # match[0]: ester oxygen, must be "O"
            if mol.GetAtomWithIdx(match[0]).GetSymbol() != "O":
                return False
            # match[1]: carbonyl carbon, must be C.
            if mol.GetAtomWithIdx(match[1]).GetSymbol() != "C":
                return False
            # For atoms indices 2..10, expect carbon.
            for idx in range(2, 11):
                if mol.GetAtomWithIdx(match[idx]).GetSymbol() != "C":
                    return False
            # Check connectivity along the chain:
            # O (at index0) should connect to C (index1)
            if not mol.GetBondBetweenAtoms(match[0], match[1]):
                return False
            # For index1, expect connection to index0 and index2.
            if not (mol.GetBondBetweenAtoms(match[1], match[2]) and mol.GetBondBetweenAtoms(match[1], match[0])):
                return False
            # For chain atoms (indices 2 to 9): each must connect to previous and next.
            for i in range(2, 10):
                if not (mol.GetBondBetweenAtoms(match[i], match[i-1]) and mol.GetBondBetweenAtoms(match[i], match[i+1])):
                    return False
            # Terminal atom (index 10): should connect only to index 9.
            if not mol.GetBondBetweenAtoms(match[10], match[9]):
                return False
            # Now check that the terminal CH3 passes the terminal check.
            if not terminal_atom_check(mol.GetAtomWithIdx(match[10]), match_set):
                return False
            return True

        elif pattern_type == "B":
            # Expected order for Pattern B:
            # match[0]: terminal CH3, should be C.
            if mol.GetAtomWithIdx(match[0]).GetSymbol() != "C":
                return False
            # match[10]: ester oxygen, must be O.
            if mol.GetAtomWithIdx(match[10]).GetSymbol() != "O":
                return False
            # Atoms indices 1 to 9 should be carbon.
            for idx in range(1, 10):
                if mol.GetAtomWithIdx(match[idx]).GetSymbol() != "C":
                    return False
            # Check connectivity along the chain:
            # Terminal at index 0: should be bonded only to index 1.
            if not mol.GetBondBetweenAtoms(match[0], match[1]):
                return False
            # For indices 1 to 8: each must bond to previous and next.
            for i in range(1, 9):
                if not (mol.GetBondBetweenAtoms(match[i], match[i-1]) and mol.GetBondBetweenAtoms(match[i], match[i+1])):
                    return False
            # Check carbonyl carbon at index 9: should be bonded to index 8 and to index 10.
            if not (mol.GetBondBetweenAtoms(match[9], match[8]) and mol.GetBondBetweenAtoms(match[9], match[10])):
                return False
            # Now check that the terminal CH3 (index 0) passes the terminal check.
            if not terminal_atom_check(mol.GetAtomWithIdx(match[0]), set(match)):
                return False
            return True
        
        else:
            return False

    # Check candidate matches for Pattern A.
    matches = mol.GetSubstructMatches(patternA)
    for match in matches:
        if check_candidate(match, "A"):
            return True, "Contains a decanoate ester functional group (matched via pattern A)"

    # Check candidate matches for Pattern B.
    matches = mol.GetSubstructMatches(patternB)
    for match in matches:
        if check_candidate(match, "B"):
            return True, "Contains a decanoate ester functional group (matched via pattern B)"

    return False, "No valid decanoate ester group found in the molecule"

# Example usage (when run as a script):
if __name__ == '__main__':
    test_examples = [
        ("O(C(CCCCCCCCC)=O)CC", "ethyl decanoate"),
        ("O(C(CCCCCCCCC)=O)CC(C)C", "isobutyl decanoate"),
        ("O=C(OC)CCCCCCCCC", "methyl decanoate"),
        # An example with a longer chain (should not be decanoate)
        ("CCCCCCCCCCCC(=O)OC", "Undecyl ester (not decanoate)"),
        # Decanoate-containing molecules from the reference list.
        ("O(C[C@@H](COC(CCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O)P(=O)(O)O", "1-decanoyl-2-stearoyl-sn-phosphatidic acid"),
    ]
    
    for smi, name in test_examples:
        result, reason = is_decanoate_ester(smi)
        print(f"Test Molecule: {name}")
        print(f"SMILES: {smi}")
        print(f"Result: {result}")
        print(f"Reason: {reason}\n")