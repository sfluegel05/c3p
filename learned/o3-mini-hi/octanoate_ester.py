"""
Classifies: CHEBI:87657 octanoate ester
"""
#!/usr/bin/env python
"""
Classifies: Octanoate Ester
Definition: Any fatty acid ester in which the carboxylic acid component is 
octanoic acid (caprylic acid). That is, the ester must contain an acyl group
with the structure: CH3–CH2–CH2–CH2–CH2–CH2–CH2–C(=O)O.
This implementation uses an explicit SMARTS pattern (with mapping numbers)
to identify the correct linear & unbranched acyl chain.
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    
    An octanoate ester (derived from octanoic acid/caprylic acid) must have its
    acyl group exactly as:
      CH3–CH2–CH2–CH2–CH2–CH2–CH2–C(=O)O
    In this SMARTS pattern, the first seven carbons represent the unbranched alkyl chain 
    (the carboxylic acid part, excluding the carbonyl carbon) and the eighth atom is 
    the carbonyl carbon.
    
    The strategy is as follows:
      1. Find all ester groups in the molecule using the general SMARTS "C(=O)O".
      2. Find all substructures matching the explicit octanoate ester SMARTS pattern:
         "[CH3:1][CH2:2][CH2:3][CH2:4][CH2:5][CH2:6][CH2:7][C:8](=O)[O]"
         (Here the mapping tag ":8" marks the carbonyl carbon of the octanoate unit.)
      3. For every ester group encountered, check if its carbonyl carbon (returned
         as the first atom in the "C(=O)O" match) is contained as the mapped carbonyl
         position (:8) in at least one octanoate match.
      4. If every ester group in the molecule meets this requirement we return True;
         otherwise, we return False together with a reason.
    
    Args:
         smiles (str): SMILES string of the molecule.
    
    Returns:
         bool: True if every ester group is derived from octanoic acid, False otherwise.
         str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General ester pattern: matches a carbonyl directly attached to an oxygen.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern, uniquify=True)
    if not ester_matches:
        return False, "No ester group found"
    
    # Define an explicit SMARTS for an octanoate group.
    # This pattern represents: CH3–CH2–CH2–CH2–CH2–CH2–CH2–C(=O)O,
    # where the carbonyl carbon is tagged as :8.
    octanoate_pattern = Chem.MolFromSmarts("[CH3:1][CH2:2][CH2:3][CH2:4][CH2:5][CH2:6][CH2:7][C:8](=O)[O]")
    octanoate_matches = mol.GetSubstructMatches(octanoate_pattern, uniquify=True)
    
    # Build a set of carbonyl atom indices that are part of an exact octanoate group.
    valid_carbonyls = set()
    for match in octanoate_matches:
        # The pattern has 8 atoms; the atom with mapping :8 is at index 7.
        if len(match) >= 8:
            valid_carbonyls.add(match[7])
    
    # Now verify every ester group in the molecule has its carbonyl carbon in valid_carbonyls.
    for match in ester_matches:
        # In the "C(=O)O" pattern, we assume the first atom is the carbonyl carbon.
        carbonyl_idx = match[0]
        if carbonyl_idx not in valid_carbonyls:
            return False, "Found an ester group whose acyl chain is not derived from octanoic acid"

    return True, "All ester groups are derived from octanoic acid (octanoate ester)."

# Example usage (testing several SMILES strings)
if __name__ == '__main__':
    test_smiles = [
        "CCCCCCCC(=O)OC[C@H](O)CO",  # 3-octanoyl-sn-glycerol: should classify True
        "CCCCCCCC(=O)OC",            # methyl octanoate: should classify True
        "O(C(=O)CCCCCCC)C(C)C",       # isopropyl octanoate: acyl chain is valid – True
        "O(C(=O)CCCCCCCCCCCC)C",      # dodecyl octanoate: acyl chain too long – False
        "CCCCCCCC(=O)OC[C@H](COP(O)(=O)OP(O)(O)=O)OC(=O)CCCCCCC",  # 1,2-dioctanoyl-sn-glycerol 3-diphosphate: True
        "CCCCCCCC(=O)O[C@H](CO)OC(=O)CCCCCCC",  # 1,2-dioctanoyl-sn-glycerol: True
        "CCCCCCCC(=O)OC",  # another valid octanoate ester (methyl octanoate)
    ]
    for sm in test_smiles:
        flag, reason = is_octanoate_ester(sm)
        print(f"SMILES: {sm}\nClassification: {flag}\nReason: {reason}\n")