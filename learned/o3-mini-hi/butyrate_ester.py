"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester
Definition: Any carboxylic ester where the acyl (acid) component is butyric acid.
That is, the ester fragment should have the form:
   [O;H0;!$([O-])]-C(=O)[CH2][CH2][C]
with the terminal carbon being isolated (only bonded to the preceding methylene).
Even when deuterated, the terminal carbon should not be connected to additional heavy atoms.
"""

from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule contains a validated butyrate ester fragment.
    The search pattern looks for an ester oxygen that is neutral and not negatively charged,
    followed by a carbonyl C, two methylene groups, and a terminal carbon.
    
    We then verify that the terminal carbon is isolated, i.e. its only heavy-atom neighbor
    (atomic number > 1) is the preceding CH2.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a validated butyrate ester fragment is found, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern for a butyrate ester fragment.
    # The pattern: an ester oxygen (not protonated or negatively charged)
    # followed by a carbonyl carbon, two methylenes, and a terminal carbon.
    # We use [C] for the terminal atom here and then check if it is an isolated methyl group.
    pattern_smarts = "[O;H0;!$([O-])]-C(=O)[CH2][CH2][C]"
    butyrate_pattern = Chem.MolFromSmarts(pattern_smarts)
    if butyrate_pattern is None:
        return False, "Invalid SMARTS pattern."
    
    # Look for all substructure matches in the molecule.
    matches = mol.GetSubstructMatches(butyrate_pattern)
    if not matches:
        return False, f"No butyrate ester fragment ({pattern_smarts}) found in the molecule."
    
    # For each match, the SMARTS match returns a tuple of indices corresponding to:
    # (ester oxygen, carbonyl carbon, first CH2, second CH2, terminal C)
    # We now post-filter to check that the terminal C is isolated.
    for match in matches:
        term_atom = mol.GetAtomWithIdx(match[4])
        # Get heavy (non-hydrogen) neighbors other than the CH2 that precedes it (match[3]).
        # Note: hydrogens are not considered heavy.
        heavy_neighbors = []
        for neighbor in term_atom.GetNeighbors():
            if neighbor.GetAtomicNum() > 1 and neighbor.GetIdx() != match[3]:
                heavy_neighbors.append(neighbor)
        # For a terminal methyl the only heavy neighbor should be the preceding CH2.
        if len(heavy_neighbors) == 0:
            return True, f"Found validated butyrate ester fragment ({pattern_smarts}) in the molecule."
    
    # If no match passed the terminality test, return False.
    return False, f"Butyrate ester fragment ({pattern_smarts}) found but not validated (terminal atom is not isolated)."

# Example usage (for testing)
if __name__ == "__main__":
    test_smiles = [
        # Examples which should be classified as butyrate esters:
        "[H][C@]12CC(=C)[C@H](CC[C@@](C)(OC(=O)CCC)[C@]([H])(O1)[C@]1([H])[C@]3([H])C[C@@]([H])(OC[C@@H]3C)[C@](C)(O)[C@]21[H])OO",  # Tritoniopsin B, rel-
        "CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCC",  # tributyrin
        "[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(OC(=O)CCC)C(=O)CO",  # cortisol 17-butyrate
        "CCCC(=O)O[C@H]1[C@H](C)O[C@H](C[C@@]1(C)O)O[C@@H]1[C@@H](C)O[C@@H](O[C@H]2[C@@H](CC=O)C[C@@H](C)[C@@H](O)\\C=C\\C=C\\C[C@@H](C)OC(=O)C[C@@H](O)[C@@H]2OC)[C@H](O)[C@H]1N(C)C",  # Leucomycin A5
        # Examples expected to be false:
        "C(CCCOC(C(CC)C)=O)CC",  # hexyl 2-methylbutanoate
        "O[C@@](C([N+](CC(=O)CCC)(C)C)([2H])[2H])(CC([O-])=O)[2H]",  # butyryl-L-carnitine-d3 test case
    ]
    
    for smi in test_smiles:
        result, reason = is_butyrate_ester(smi)
        print(f"SMILES: {smi}\nResult: {result} | Reason: {reason}\n")