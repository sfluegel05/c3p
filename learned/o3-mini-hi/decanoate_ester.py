"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: decanoate ester
Definition:
 A decanoate ester is defined as a fatty acid ester resulting from the formal condensation 
 of decanoic acid [CH3(CH2)8COOH] with an alcohol or phenol. In the acyl group (decanoate), 
 the carbonyl carbon (C=O) plus the alkyl chain must consist of exactly 10 carbons 
 (i.e. CH3-(CH2)8-C(=O)-).
 
Our improved approach:
 1. Define two SMARTS patterns for the acyl decanoate fragment:
    Pattern A: "OC(=O)CCCCCCCCC" (ester oxygen is on the left).
    Pattern B: "CCCCCCCCC(=O)O" (ester oxygen is on the right).
 2. For each match (which should return 11 atoms in the match), check that the terminal 
    alkyl carbon (index 10 in Pattern A, index 0 in Pattern B) is a terminal CH3:
      • It must have only one bonded carbon (its neighbor in the chain).
      • Its total number of attached hydrogens (implicit+explicit) should be 3.
 3. If any candidate match passes these extra checks then we classify the molecule as 
    containing a decanoate ester.
 4. Otherwise we return that no valid decanoate ester group was found.
"""

from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester must contain an ester linkage where the decanoate acyl group 
    (CH3(CH2)8CO–; exactly 10 carbons) is present.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a valid decanoate ester group, False otherwise.
        str : Reason for classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns:
    # Pattern A: ester oxygen on the left: O-C(=O)-CCCCCCCCC (match should have 11 atoms)
    patternA = Chem.MolFromSmarts("OC(=O)CCCCCCCCC")
    # Pattern B: ester oxygen on the right: CCCCCCCCC(=O)O (match should also have 11 atoms)
    patternB = Chem.MolFromSmarts("CCCCCCCCC(=O)O")
    
    # Helper function to check that a candidate terminal alkyl carbon (expected to be a CH3)
    # is indeed terminal and unbranched.
    def terminal_check(atom, match_indices):
        """
        Ensures that the given atom is a terminal methyl carbon.
        It must have exactly one carbon neighbor and exactly three hydrogens.
        """
        # First, it must be a carbon.
        if atom.GetSymbol() != "C":
            return False
        # The degree (number of bonds excluding hydrogens) should be 1 for a CH3 terminal.
        if atom.GetDegree() != 1:
            return False
        # Additionally, check the total number of hydrogens (implicit+explicit).
        if atom.GetTotalNumHs() != 3:
            return False
        # Also ensure that its only bonded carbon is part of the candidate match.
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        in_match = [nbr for nbr in carbon_neighbors if nbr.GetIdx() in match_indices]
        if len(carbon_neighbors) == 1 and len(in_match) == 1:
            return True
        return False

    # Check candidate matches for Pattern A.
    matches = mol.GetSubstructMatches(patternA)
    for match in matches:
        # We expect exactly 11 atoms in a valid match.
        if len(match) != 11:
            continue
        # For Pattern A, by our convention, the terminal alkyl carbon should be at index 10.
        terminal_idx = match[10]
        terminal_atom = mol.GetAtomWithIdx(terminal_idx)
        if terminal_check(terminal_atom, match):
            return True, "Contains a decanoate ester functional group (matched via pattern A)"
    
    # Now check candidate matches for Pattern B.
    matches = mol.GetSubstructMatches(patternB)
    for match in matches:
        if len(match) != 11:
            continue
        # For Pattern B, the terminal alkyl carbon is expected to be at index 0.
        terminal_idx = match[0]
        terminal_atom = mol.GetAtomWithIdx(terminal_idx)
        if terminal_check(terminal_atom, match):
            return True, "Contains a decanoate ester functional group (matched via pattern B)"
    
    # If no valid decanoate ester fragment is found, return False.
    return False, "No valid decanoate ester group found in the molecule"

# Example usage:
if __name__ == '__main__':
    # Test SMILES examples (these include correct decanoate esters and other controls)
    test_smiles_list = [
        ("O(C(CCCCCCCCC)=O)CC", "ethyl decanoate"),
        ("O(C(CCCCCCCCC)=O)CC(C)C", "isobutyl decanoate"),
        ("O=C(OC)CCCCCCCCC", "methyl decanoate"),
        # This one is not decanoate (longer alkyl chain, for example)
        ("CCCCCCCCCCCC(=O)OC", "Undecyl ester (not decanoate)"),
        # One of the reference compounds:
        ("O(C[C@@H](COC(CCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCCC)=O)P(=O)(O)O", "1-decanoyl-2-stearoyl-sn-phosphatidic acid"),
    ]
    
    for smi, name in test_smiles_list:
        result, reason = is_decanoate_ester(smi)
        print(f"Test Molecule: {name}")
        print(f"SMILES: {smi}")
        print(f"Result: {result}")
        print(f"Reason: {reason}\n")