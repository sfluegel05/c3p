"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside

Definition: Any D-galactoside having beta-configuration at its anomeric centre.
A beta-D-galactoside here is defined as a molecule that contains a complete beta-D-galactopyranoside fragment.
A complete beta-D-galactopyranoside fragment is defined as a pyranose ring (six‐membered ring with one oxygen
and five carbons) in which an external substituent (the aglycone) is attached via an oxygen to the anomeric carbon;
the anomeric carbon is drawn with beta configuration ([C@@H]), and the ring atoms have the expected relative arrangement.
In beta-D-galactopyranose the exocyclic CH2OH is attached at C5.
Because many alternative representations exist, this approach is strict and may fail if stereo is not defined.
"""

from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    
    The function first parses the SMILES and assigns stereochemistry. It uses a refined SMARTS pattern
    that looks for an external substituent attached via an oxygen to an anomeric carbon in beta configuration.
    The pattern enforces a complete six-membered pyranose ring (with one ring oxygen and five carbons)
    with the exocyclic CH2OH group bonded at the appropriate ring carbon.
    
    The SMARTS pattern we use here is:
      [*]O[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H](CO)O1
    This denotes:
      [*]O       : any aglycone attached via an oxygen
      [C@@H]1   : the anomeric carbon in beta configuration (start ring number 1)
         O      : the ring oxygen
         [C@H](O): C2 with an OH
         [C@@H](O): C3 with an OH
         [C@H](O): C4 with an OH
         [C@H](CO): C5 with an exocyclic CH2OH (note the CH2OH group, not in the ring)
      O1         : closes the ring (the ring oxygen already included)
    
    In addition to matching the SMARTS, we verify that the atoms in the sugar ring (the atoms with indices
    corresponding to the anomeric carbon and the next four atoms plus the oxygen closing the ring)
    actually form a ring of size 6.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a complete beta-D-galactoside fragment, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is computed (this is important for chiral centers in sugars).
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define a refined SMARTS pattern for a beta-D-galactopyranoside.
    # This pattern requires an external substituent (any atom, [*]) attached via an oxygen
    # to a beta-configured anomeric carbon which is part of a six-membered ring with exactly one oxygen.
    smarts_beta_gal = "[*]O[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H](CO)O1"
    patt = Chem.MolFromSmarts(smarts_beta_gal)
    if patt is None:
        return None, "Error creating SMARTS pattern"
    
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return False, "Beta-D-galactoside moiety not found (no SMARTS match)"
    
    # Get ring information for further filtering.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # For each match, check that the sugar ring is complete.
    # In our SMARTS, the ring atoms are in positions 1 through 6 of the match.
    # (Match tuple: index0 = aglycone O, then index1 = anomeric carbon, index2 = ring oxygen, etc.)
    for match in matches:
        # Get the set of atom indices corresponding to the ring in the SMARTS.
        # The pattern uses ring label "1" starting at the anomeric carbon (match[1]) and closing at match[6].
        # Thus we consider atoms match[1] ... match[6].
        if len(match) < 7:
            continue  # Skip if the match is incomplete.
        ring_atoms = set(match[1:7])
        # Now check if these atoms form a ring of size 6.
        for ring in ring_info:
            if len(ring) == 6 and ring_atoms.issubset(ring):
                # We found a complete sugar ring that matches the expected pattern.
                return True, "Contains beta-D-galactoside moiety based on refined SMARTS matching and ring verification"
    
    return False, "No complete beta-D-galactoside sugar ring found despite SMARTS match"


# Example usage (for testing):
if __name__ == "__main__":
    # Test with a confirmed beta-D-galactoside: methyl beta-D-galactoside
    test_smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"
    result, reason = is_beta_D_galactoside(test_smiles)
    print(test_smiles, "->", result, "|", reason)