"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
Our revised strategy:
  • We identify a glycosidic bond by finding an exocyclic oxygen (“-O–”) attached to a chiral 
    anomeric carbon that is drawn with beta stereochemistry (using the ‘@@’ notation or chirality 
    tag without insisting on an explicit hydrogen).
  • We use three SMARTS patterns that differ only in the allowed substituent on the anomeric carbon.
    In each pattern, we require a six-membered (pyranose) sugar ring:
      Pattern A: Requires an explicit CH2OH at the exocyclic position.
      Pattern B: Allows an extra oxygen (e.g. acetylation) at that same position.
      Pattern C: More general, allowing any saturated substituent.
  • Finally, among all SMARTS matches we verify that at least one of them “contains” an entire six‐membered ring.
  
Note: As before, we note that there exist many valid SMILES representations of sugars and glycosides.
While our SMARTS patterns aim to include many common cases for beta-D-glucosides, some molecules may be mis‐classified.
"""

from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    Our revised approach is to search for an O–glycosidic linkage to a sugar ring whose anomeric
    carbon is drawn with beta stereochemistry (using the '@@' chirality descriptor).
    
    We use three SMARTS patterns that match an external oxygen (which must carry at least one H)
    bound to a chiral anomeric carbon (expressed as [C@@] so that we do not force an explicit hydrogen)
    that is part of a six-membered ring similar to glucopyranose.
    
    Patterns:
      Pattern A: Expects a CH2OH substituent at the anomeric carbon.
                  "*-O[C@@]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
      Pattern B: Allows for additional oxygenation at that substituent (e.g. an acylated group).
                  "*-O[C@@]1O[C@H](C(O*))[C@@H](O)[C@H](O)[C@H]1O"
      Pattern C: More general: any saturated (sp3) substituent is allowed.
                  "*-O[C@@]1O[C@H]([CX4])[C@@H](O)[C@H](O)[C@H]1O"
                  
    After substructure matching, we require that one of the matches “covers” a complete six‐membered ring,
    based on the molecule’s ring information.
    
    Args:
        smiles (str): SMILES string.
    Returns:
        bool: True if the molecule is classified as a beta-D-glucoside, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns.
    # The leading "*" permits any substituent at the glycosidic donor.
    patternA_smarts = "*-O[C@@]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    patternB_smarts = "*-O[C@@]1O[C@H](C(O*))[C@@H](O)[C@H](O)[C@H]1O"
    patternC_smarts = "*-O[C@@]1O[C@H]([CX4])[C@@H](O)[C@H](O)[C@H]1O"
    
    patA = Chem.MolFromSmarts(patternA_smarts)
    patB = Chem.MolFromSmarts(patternB_smarts)
    patC = Chem.MolFromSmarts(patternC_smarts)
    
    if not (patA and patB and patC):
        return False, "Error in SMARTS pattern definition"
    
    # Gather all matches (each match is a tuple of atom indices in mol)
    matches = []
    matches.extend(mol.GetSubstructMatches(patA))
    matches.extend(mol.GetSubstructMatches(patB))
    matches.extend(mol.GetSubstructMatches(patC))
    
    if not matches:
        return False, "No beta-D-glucoside substructure (with beta anomeric centre) found."
    
    # Now verify that at least one of these matches covers a full six-membered ring (pyranose ring).
    # We use the ring information from the molecule (each ring is a tuple of atom indices).
    ring_info = mol.GetRingInfo().AtomRings()
    valid_match_found = False
    # Remove duplicate matches by using a set of tuples.
    unique_matches = {tuple(match) for match in matches}
    for match in unique_matches:
        # For each match, check each ring: if a ring of size 6 is completely contained in the match,
        # assume that the match represents a full sugar ring.
        for ring in ring_info:
            if len(ring) == 6:
                if set(ring).issubset(set(match)):
                    valid_match_found = True
                    break
        if valid_match_found:
            break
    
    if valid_match_found:
        return True, "Molecule contains a beta-D-glucoside moiety with beta-anomeric configuration."
    else:
        return False, "Matching fragment was found but does not cover a complete six-membered sugar ring."

# (Optional) Example tests:
if __name__ == '__main__':
    # A few examples are provided; note that the full test set is extensive.
    test_examples = [
        # True positives (known beta-D-glucosides)
        "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC",  # beta-D-glucosyl-N-eicosanoylsphingosine
        "[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C",      # oleuropein
        "OC[C@H]1O[C@@H](OCCc2ccc(O)c(O)c2)[C@H](OC(=O)Cc2ccc(O)cc2)[C@@H](O)[C@@H]1O",  # ternstroside B
        # False negatives (should be matched but previously missed) examples:
        "[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)COC(=O)C)O)O)O)O[C@@]23[C@]4([C@@](C2)([C@@]5(O[C@]4(O[C@]3(C5)C)[H])O)[H])COC(C6=CC=CC=C6)=O",  # 6'-O-acetylpaeoniflorin (expect True)
    ]
    for s in test_examples:
        res, reason = is_beta_D_glucoside(s)
        print(f"SMILES: {s}\nResult: {res}, Reason: {reason}\n")