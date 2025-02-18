"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
This code defines several SMARTS patterns – including one for N‐acetylated sugars – that force an explicit hydrogen on the anomeric carbon (using [C@@H])
to bias the match toward beta stereochemistry. Then it verifies that the matching fragment covers a complete six‐membered (pyranose) ring.
Note: Some valid representations may still be missed and some false matches occur because of the many ways sugar moieties are drawn.
"""

from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES.
    Approach:
      • Define several SMARTS patterns that capture an O-glycosidic linkage from the glycosidic donor (external O)
        to a sugar ring whose anomeric carbon is drawn with beta stereochemistry.
      • Patterns A–C cover a range of substituents at the anomeric carbon while forcing the chiral centre
        to be written as [C@@H] (explicit hydrogen) so that the beta configuration is enforced.
      • Pattern D is added to capture the common N-acetyl case (e.g. in beta-D-GlcpNAc).
      • Afterwards, all substructure matches are checked to ensure that one of them contains an entire six-membered ring.
    
    Args:
        smiles (str): input SMILES string.
    
    Returns:
        (bool, str): True with a success message if a beta-D-glucoside moiety is detected; else False and a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define several SMARTS patterns.
    # The leading "*" allows any substituent on the donor side.
    # We use [C@@H] (rather than [C@@]) so that an explicit hydrogen is required at the anomeric carbon.
    patternA_smarts = "*-O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    # Pattern B allows an alternative substituent at the anomeric carbon; for example if it is oxygenated (e.g. acetylated)
    patternB_smarts = "*-O[C@@H]1O[C@H](C(O*))[C@@H](O)[C@H](O)[C@H]1O"
    # Pattern C is more general: any saturated substituent at the anomeric position.
    patternC_smarts = "*-O[C@@H]1O[C@H]([CX4])[C@@H](O)[C@H](O)[C@H]1O"
    # Pattern D targets N-acetyl substitution on the sugar (e.g. GlcpNAc); notice the NC(=O)C group.
    patternD_smarts = "*-O[C@@H]1O[C@H](NC(=O)C)[C@@H](O)[C@H](O)[C@H]1O"
    
    patterns_smarts = [patternA_smarts, patternB_smarts, patternC_smarts, patternD_smarts]
    
    patterns = []
    for smarts in patterns_smarts:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            return False, "Error in SMARTS pattern definition"
        patterns.append(pat)
    
    # Gather all unique matches (each match is a tuple of atom indices)
    matches = set()
    for pat in patterns:
        for m in mol.GetSubstructMatches(pat):
            matches.add(m)
    
    if not matches:
        return False, "No beta-D-glucoside substructure (with beta anomeric centre) found."
    
    # Use the ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    valid_match_found = False
    for match in matches:
        # For each match, check whether any full six-membered ring is contained.
        for ring in ring_info:
            if len(ring) == 6 and set(ring).issubset(set(match)):
                valid_match_found = True
                break
        if valid_match_found:
            break
    
    if valid_match_found:
        return True, "Molecule contains a beta-D-glucoside moiety with beta-anomeric configuration."
    else:
        return False, "Fragment matching beta-D-glucoside found, but no complete six-membered sugar ring was covered."

# (Optional) If run as main, test some examples.
if __name__ == '__main__':
    # A few test examples are provided.
    test_examples = [
        # True positives (known beta-D-glucosides)
        ("CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC",
         "beta-D-glucosyl-N-eicosanoylsphingosine"),
        ("[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C",
         "oleuropein"),
        ("OC[C@H]1O[C@@H](OCCc2ccc(O)c(O)c2)[C@H](OC(=O)Cc2ccc(O)cc2)[C@@H](O)[C@@H]1O",
         "ternstroside B"),
        # One example that was previously missed (beta-D-GlcpNAc-(1->6)-GalNAc-ol)
        ("O1[C@@H](OC[C@@H](O)[C@H](O)[C@H](O)[C@@H](NC(=O)C)CO)[C@H](NC(=O)C)[C@@H](O)[C@H](O)[C@H]1CO",
         "beta-D-GlcpNAc-(1->6)-GalNAc-ol")
    ]
    
    for s, name in test_examples:
        res, reason = is_beta_D_glucoside(s)
        print(f"Name: {name}\nSMILES: {s}\nResult: {res}, Reason: {reason}\n")