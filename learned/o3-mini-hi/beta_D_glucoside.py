"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
Our strategy:
  • We search for a beta-D-glucopyranoside fragment. In our approach the beta-configuration
    is enforced by requiring that the glycosidic oxygen is bound to an anomeric carbon drawn with 
    the beta stereochemistry ([C@@H]). In a typical fragment the sugar ring is a pyranose (six-membered)
    ring with one ring oxygen.
  • To allow for chemical variability (acylation, substitutions on the CH2OH group, etc.) we use three 
    slightly different SMARTS patterns that only differ in the substituent at the anomeric carbon.
      Pattern A: a simple CH2OH substituent.
      Pattern B: allows for additional oxygenation (e.g. acetylation) on that carbon.
      Pattern C: more general (any sp3 carbon is allowed).
  • Finally, we require that the atoms in any substructure match include a six-membered ring,
    as determined by the molecule’s ring information.
    
Note: Because many representations of sugars and glycosidic bonds exist in SMILES, this approach
may mis‐classify some molecules. Our patterns were chosen in an attempt to capture as much of the beta‐D‐glucoside
class as possible while avoiding false positives.
"""

from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    Our definition is that the molecule should contain a glucopyranose (six-membered sugar) 
    unit with an exocyclic oxygen substituent attached to an anomeric carbon that is drawn 
    with beta stereochemistry ([C@@H] in the SMILES). 
    
    We use three SMARTS patterns that differ only in what is allowed at the position 
    immediately bonded to the anomeric carbon. The patterns are:
      Pattern A (simple CH2OH):  
           [!H0]-O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O
      Pattern B (allows extra oxygenation e.g. acylation at the CH2OH):  
           [!H0]-O[C@@H]1O[C@H](C(O*))[C@@H](O)[C@H](O)[C@H]1O
      Pattern C (more general; any sp3-carbon at that position):  
           [!H0]-O[C@@H]1O[C@H]([CX4])[C@@H](O)[C@H](O)[C@H]1O

    In all cases, [!H0] ensures that the oxygen substituent actually bears at least one hydrogen
    (i.e. it is not a deprotonated oxygen), and the [C@@H] in the fragment forces beta stereochemistry at the anomeric centre.
    
    Additionally, we check that at least one substructure match “covers” a six-membered ring.
    
    Args:
        smiles (str): SMILES string.
        
    Returns:
        bool: True if the molecule is classified as a beta-D-glucoside, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define three SMARTS patterns; note that extra substituents on the CH2OH group are allowed.
    patternA_smarts = "[!H0]-O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    patternB_smarts = "[!H0]-O[C@@H]1O[C@H](C(O*))[C@@H](O)[C@H](O)[C@H]1O"
    patternC_smarts = "[!H0]-O[C@@H]1O[C@H]([CX4])[C@@H](O)[C@H](O)[C@H]1O"
    
    patA = Chem.MolFromSmarts(patternA_smarts)
    patB = Chem.MolFromSmarts(patternB_smarts)
    patC = Chem.MolFromSmarts(patternC_smarts)
    
    if not (patA and patB and patC):
        return False, "Error in SMARTS pattern definition"
    
    # Find all matches from any of the three patterns.
    matches = []
    matches.extend(mol.GetSubstructMatches(patA))
    matches.extend(mol.GetSubstructMatches(patB))
    matches.extend(mol.GetSubstructMatches(patC))
    
    if not matches:
        return False, "No beta-D-glucoside substructure (with beta anomeric centre) found."
    
    # Verify that at least one match corresponds to a six-membered (pyranose) ring.
    ring_info = mol.GetRingInfo().AtomRings()
    valid_match_found = False
    # Remove duplicate matches by turning them into tuples.
    unique_matches = {tuple(match) for match in matches}
    for match in unique_matches:
        # For each match, check every ring in the molecule.
        for ring in ring_info:
            if len(ring) == 6:
                # If the entire ring is a subset of the matched atoms, we assume that the sugar ring was found.
                if set(ring).issubset(set(match)):
                    valid_match_found = True
                    break
        if valid_match_found:
            break
    
    if valid_match_found:
        return True, "Molecule contains a beta-D-glucoside moiety with beta-anomeric configuration."
    else:
        return False, "Matching fragment not part of a complete six-membered sugar ring."

# (Optional) Example tests:
if __name__ == '__main__':
    test_examples = [
        # True positives (examples known to belong to the beta-D-glucoside class)
        "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC",  # beta-D-glucosyl-N-eicosanoylsphingosine
        "[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C",  # oleuropein
        "OC[C@H]1O[C@@H](OCCc2ccc(O)c(O)c2)[C@H](OC(=O)Cc2ccc(O)cc2)[C@@H](O)[C@@H]1O",  # ternstroside B
        # (Additional examples from the list could be added here)
    ]
    for s in test_examples:
        res, reason = is_beta_D_glucoside(s)
        print(f"SMILES: {s}\nResult: {res}, Reason: {reason}\n")