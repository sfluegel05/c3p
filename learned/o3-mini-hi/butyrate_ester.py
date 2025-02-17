"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester
Definition: Any carboxylic ester where the acid (acyl) component is butyric acid.
That is the acyl group is butyrate (CH3CH2CH2CO–); we look for an ester fragment:
    –O–C(=O)–CH2–CH2–C
and then verify that the last carbon is terminal (i.e. its only neighbor is the previous methylene).
This improves on a naïve SMARTS by:
  • Allowing for deuterated terminal methyl groups.
  • Checking that the putative butyrate acyl portion is not “embedded” in a larger system.
"""

from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule contains a butyrate ester fragment (i.e. one in which the acid part is butyric acid).
    Here we search for a substructure that looks like:
       [O;H0;!$([O-])]-C(=O)[CH2][CH2][C]
    and we then verify that the terminal C is in fact terminal (has only the expected connection).
    
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
    
    # Define a SMARTS pattern for a butyrate ester fragment.
    # Instead of insisting that the terminal group is exactly [CH3;D1],
    # we require a sequence of atoms: ester oxygen bonded to a carbonyl carbon,
    # then two methylenes, and finally a carbon (which may be CH3 or deuterated).
    # Later we check that the last carbon is terminal.
    pattern_smarts = "[O;H0;!$([O-])]-C(=O)[CH2][CH2][C]"
    butyrate_pattern = Chem.MolFromSmarts(pattern_smarts)
    if butyrate_pattern is None:
        return False, "Invalid SMARTS pattern."
    
    # Look for all substructure matches in the molecule.
    # The match tuple will have indices corresponding to:
    # (ester O, carbonyl C, first CH2, second CH2, terminal C)
    matches = mol.GetSubstructMatches(butyrate_pattern)
    if not matches:
        return False, "No butyrate ester fragment (%s) found in the molecule." % pattern_smarts

    # Post-filter each match to verify that the terminal carbon is really terminal.
    # (This check helps catch cases where the pattern might be hit in contorted settings,
    # or where the terminal atom is not acting as a methyl group from a free butyric acid unit.)
    for match in matches:
        # Get the atom corresponding to the terminal carbon
        term_atom = mol.GetAtomWithIdx(match[4])
        # Use GetTotalDegree() to count all bonds (explicit+implicit). For a terminal methyl,
        # its only connection should be to the previous methylene.
        if term_atom.GetTotalDegree() == 1:
            return True, "Found validated butyrate ester fragment (%s) in the molecule." % pattern_smarts

    # If none of the matches pass the terminality check, then reject.
    return False, "Butyrate ester fragment (%s) found but not validated (terminal atom not isolated)." % pattern_smarts

# Example usage (for testing)
if __name__ == "__main__":
    test_smiles = [
        # True positives – these molecules contain the expected butyrate ester fragment.
        "CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCC",   # tributyrin
        "CCCC(=O)OCC",                        # ethyl butyrate
        "[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(OC(=O)CCC)C(=O)CO", # cortisol 17-butyrate
        # False negatives (should not classify as butyrate ester):
        "C(CCCOC(C(CC)C)=O)CC",                # hexyl 2-methylbutanoate (not butyric acid)
        "O[C@@](C([N+](CC(=O)CCC)(C)C)([2H])[2H])(CC([O-])=O)[2H]",  # butyryl-L-carnitine-d3 test case
    ]
    
    for smi in test_smiles:
        result, reason = is_butyrate_ester(smi)
        print(f"SMILES: {smi}\nResult: {result} | Reason: {reason}\n")