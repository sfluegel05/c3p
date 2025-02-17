"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: Glycosphingolipid
A glycosphingolipid is defined as a glycolipid that is a carbohydrate‐containing derivative 
of a sphingoid or ceramide, with the carbohydrate residue attached by a glycosidic linkage 
to O-1 of the sphingoid.
"""

from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    The strategy is to (1) confirm that a carbohydrate (sugar) ring is present,
    (2) detect a sphingoid/ceramide motif (here, a pattern corresponding to an amide 
        attached to a chiral carbon with a hydroxymethyl branch, common in sphingosine),
    and (3) verify that the sugar is attached to the sphingoid via an oxygen (glycosidic linkage).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for the classification
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Look for a carbohydrate “sugar” moiety.
    # Here we use a SMARTS pattern for a typical six-membered (pyranose) sugar.
    sugar_smarts = "O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)C1"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No carbohydrate (sugar) moiety detected"
    
    # Step 2: Look for a sphingoid/ceramide backbone.
    # A common motif in sphingolipids is an amide-linked chain containing a chiral center with a 
    # hydroxymethyl branch. In many glycosphingolipids the sphingoid base shows up with a fragment 
    # like "N[C@@H](CO...". We use this as a heuristic.
    sphingo_smarts = "N[C@@H](CO)"
    sphingo_pattern = Chem.MolFromSmarts(sphingo_smarts)
    sphingo_matches = mol.GetSubstructMatches(sphingo_pattern)
    if not sphingo_matches:
        return False, "No sphingoid/ceramide backbone detected"
    
    # Step 3 (optional but useful): Check for a glycosidic linkage.
    # We expect that the sugar moiety is attached by an oxygen atom to the sphingoid base.
    # We collect all atom indices in the sugar matches and in the sphingoid matches.
    sugar_atoms = set()
    for match in sugar_matches:
        sugar_atoms.update(match)
    sphingo_atoms = set()
    for match in sphingo_matches:
        sphingo_atoms.update(match)
    
    glycosidic_link_found = False
    # Iterate over all bonds in the molecule. Look for a bond that connects a sugar atom to a sphingoid atom.
    # To increase our confidence that the connection is via an oxygen,
    # we require that at least one end of the bond is an oxygen atom.
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        if a.GetAtomicNum() == 8 or b.GetAtomicNum() == 8:
            if ((a.GetIdx() in sugar_atoms and b.GetIdx() in sphingo_atoms) or 
                (b.GetIdx() in sugar_atoms and a.GetIdx() in sphingo_atoms)):
                glycosidic_link_found = True
                break
                
    if not glycosidic_link_found:
        return False, "Sugar and sphingoid moieties are not glycosidically linked"
    
    return True, "Molecule contains a glycosphingolipid structure: carbohydrate moiety glycosidically linked to a sphingoid/ceramide backbone"

# Example usage (uncomment to test):
# smiles_example = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_glycosphingolipid(smiles_example)
# print(result, reason)