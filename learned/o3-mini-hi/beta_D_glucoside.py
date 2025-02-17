"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
Our strategy:
  • Look for a beta-D-glucopyranoside substructure. A typical fragment is drawn as:
        [!H0]-O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O
    which requires that the glycosidic bond is to an anomeric carbon drawn with the beta (i.e. @@) configuration.
  • Because substitutions on the CH2OH (i.e. "CO") unit are common, we also allow a variant
        [!H0]-O[C@@H]1O[C@H](C(O*))[C@@H](O)[C@H](O)[C@H]1O
  • Finally, to be more flexible the third pattern replaces the explicit “CO” with a general sp3-carbon ([CX4]),
        yielding: [!H0]-O[C@@H]1O[C@H]([CX4])[C@@H](O)[C@H](O)[C@H]1O
  • In addition we check that the atoms forming the sugar fragment belong to a 6-membered ring (i.e. a pyranose).
"""

from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside (in our definition) is a molecule that contains a glucopyranose
    (six-membered sugar) unit attached via an oxygen to an aglycone and in which
    the anomeric carbon is drawn with beta stereochemistry.
    
    We use three SMARTS patterns:
      Pattern A (classic unmodified beta-D-glucopyranoside):
          [!H0]-O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O
      Pattern B (allows substitution on the CH2OH group):
          [!H0]-O[C@@H]1O[C@H](C(O*))[C@@H](O)[C@H](O)[C@H]1O
      Pattern C (more general, allowing any sp3 carbon at that position):
          [!H0]-O[C@@H]1O[C@H]([CX4])[C@@H](O)[C@H](O)[C@H]1O

    Additionally, we only accept a match if the atoms belonging to the sugar fragment
    form a six-membered ring.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule is classified as a beta-D-glucoside, False otherwise.
        str : Reason for classification.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    patternA_smarts = "[!H0]-O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    patternB_smarts = "[!H0]-O[C@@H]1O[C@H](C(O*))[C@@H](O)[C@H](O)[C@H]1O"
    patternC_smarts = "[!H0]-O[C@@H]1O[C@H]([CX4])[C@@H](O)[C@H](O)[C@H]1O"
    
    patternA = Chem.MolFromSmarts(patternA_smarts)
    patternB = Chem.MolFromSmarts(patternB_smarts)
    patternC = Chem.MolFromSmarts(patternC_smarts)
    
    if not (patternA and patternB and patternC):
        return False, "Error in SMARTS pattern definition"
    
    # Collect all substructure matches from any pattern
    matches = []
    matches.extend(mol.GetSubstructMatches(patternA))
    matches.extend(mol.GetSubstructMatches(patternB))
    matches.extend(mol.GetSubstructMatches(patternC))
    
    if not matches:
        return False, "No beta-D-glucoside substructure (with beta anomeric centre) found."
    
    # For each match, verify that the atoms of the sugar fragment form a six-membered ring.
    # (A typical beta-D-glucopyranoside has a pyranose ring – a ring of 6 atoms.)
    ring_info = mol.GetRingInfo().AtomRings()
    valid_match_found = False
    for match in set(tuple(m) for m in matches):  # use set to remove duplicate matches
        # Check every ring in the molecule to see if any are a 6-membered ring contained in the match.
        for ring in ring_info:
            if len(ring) == 6 and set(ring).issubset(set(match)):
                valid_match_found = True
                break
        if valid_match_found:
            break

    if valid_match_found:
        return True, "Molecule contains a beta-D-glucoside moiety with beta-anomeric configuration."
    else:
        return False, "No valid beta-D-glucoside (6-membered glucopyranose ring with beta configuration) found."

# Example tests (you can uncomment these lines to run some tests)
# test_smiles = [
#     "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC",  # beta-D-glucosyl-N-eicosanoylsphingosine
#     "[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C"  # oleuropein
# ]
# for s in test_smiles:
#     res, reason = is_beta_D_glucoside(s)
#     print(s, "=>", res, reason)