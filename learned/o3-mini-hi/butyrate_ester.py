"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester
Definition: Any carboxylic ester where the acyl (acid) component is butyric acid.
That is, the ester fragment should have the form:
   [O;H0;!$([O-])]-C(=O)[CH2][CH2][C]
and further the terminal carbon (the last C) must be isolated:
it should only be bonded to the preceding methylene (i.e. be a terminal methyl group),
even when deuterated.
"""

from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule contains a validated butyrate ester fragment.
    We first search for the substructure pattern that looks like:
       [O;H0;!$([O-])]-C(=O)[CH2][CH2][C]
    This corresponds to the ester oxygen, the carbonyl carbon, two methylene groups,
    and a terminal carbon. We then verify that the terminal carbon is in fact terminal,
    i.e. it has only one explicit heavy-atom neighbor.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a validated butyrate ester fragment is found, False otherwise.
        str: Explanation of the classification.
    """
    # Attempt to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define and compile a SMARTS pattern for a butyrate ester fragment.
    # The pattern corresponds to:
    # Ester oxygen that is not protonated and not negatively charged,
    # followed by a carbonyl C, two CH2 groups, and finally a carbon.
    pattern_smarts = "[O;H0;!$([O-])]-C(=O)[CH2][CH2][C]"
    butyrate_pattern = Chem.MolFromSmarts(pattern_smarts)
    if butyrate_pattern is None:
        return False, "Invalid SMARTS pattern."
    
    # Look for all substructure matches in the molecule.
    # The match tuple indices are:
    # (ester oxygen, carbonyl C, first CH2, second CH2, terminal C).
    matches = mol.GetSubstructMatches(butyrate_pattern)
    if not matches:
        return False, f"No butyrate ester fragment ({pattern_smarts}) found in the molecule."
    
    # Post-filter each match to verify that the terminal carbon is really a terminal methyl group.
    # For an isolated CH3, its explicit heavy atom degree should be exactly 1 (i.e. it is only bonded
    # to the preceding CH2 in the butyrate fragment).
    for match in matches:
        term_atom = mol.GetAtomWithIdx(match[4])
        if term_atom.GetDegree() == 1:  # Check only explicit (heavy-atom) neighbors.
            return True, f"Found validated butyrate ester fragment ({pattern_smarts}) in the molecule."
    
    # If none of the matches pass the terminality check, then reject.
    return False, f"Butyrate ester fragment ({pattern_smarts}) found but not validated (terminal atom not isolated)."

# Example usage (for testing)
if __name__ == "__main__":
    test_smiles = [
        # Should be classified as butyrate esters.
        "[H][C@]12CC(=C)[C@H](CC[C@@](C)(OC(=O)CCC)[C@]([H])(O1)[C@]1([H])[C@]3([H])C[C@@]([H])(OC[C@@H]3C)[C@](C)(O)[C@]21[H])OO",  # Tritoniopsin B, rel-
        "CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCC",  # tributyrin
        "[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(OC(=O)CCC)C(=O)CO",  # cortisol 17-butyrate
        "CCCC(=O)O[C@H]1[C@H](C)O[C@H](C[C@@]1(C)O)O[C@@H]1[C@@H](C)O[C@@H](O[C@H]2[C@@H](CC=O)C[C@@H](C)[C@@H](O)\\C=C\\C=C\\C[C@@H](C)OC(=O)C[C@@H](O)[C@@H]2OC)[C@H](O)[C@H]1N(C)C",  # Leucomycin A5
        # Examples that are not butyrate esters (should be false).
        "C(CCCOC(C(CC)C)=O)CC",  # hexyl 2-methylbutanoate
        "O[C@@](C([N+](CC(=O)CCC)(C)C)([2H])[2H])(CC([O-])=O)[2H]",  # butyryl-L-carnitine-d3 test case
    ]
    
    for smi in test_smiles:
        result, reason = is_butyrate_ester(smi)
        print(f"SMILES: {smi}\nResult: {result} | Reason: {reason}\n")