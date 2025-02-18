"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: Parent monosaccharide.
A parent monosaccharide is defined heuristically as either:
  (a) An open‐chain polyhydroxy aldehyde (or ketone) with a terminal aldehyde (or internal ketone)
      if sufficient –OH groups are present, or
  (b) A cyclic sugar (5‐ or 6‐membered ring containing oxygen) with several –OH groups.
Only molecules consisting solely of carbon and oxygen (with formal charge 0 on all atoms)
and having 3–9 carbon atoms and a single fragment are accepted.
If no typical sugar motifs are found, the compound is rejected.
Note: This is a heuristic; if the task is too difficult we may return (None, None).
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a parent monosaccharide based on its SMILES string.
    Improved heuristics:
      1. Only allow molecules with C and O atoms (H is implicit).
      2. Must be a single fragment and all atoms must have formal charge zero.
      3. Carbon atoms count must be between 3 and 9.
      4. Count hydroxyl groups using the SMARTS "[OX2H]".
      5. Check for a typical cyclic sugar ring: a 5- or 6-membered ring containing at least one oxygen.
         – If found, require at least 3 hydroxyl groups and an OH/C ratio ≥ 0.5.
      6. If no ring is found, seek an open‐chain carbonyl:
         – For an aldose: a terminal aldehyde pattern "[H]C(=O)" is required and expect roughly (C-1) –OH groups.
         – For a ketose: require an internal carbonyl pattern "C(=O)C" and expect roughly (C-2) –OH groups.
      7. Reject if any extraneous elements (other than C and O) are present.
      8. Reject molecules with any non‐zero formal charge.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a parent monosaccharide, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "SMILES parsing failed"

    # 1. Check connectivity: single fragment only.
    frags = Chem.GetMolFrags(mol, asMols=False)
    if len(frags) > 1:
        return False, "Multiple fragments detected; not a single monosaccharide unit"

    # 2. Check that every atom is either carbon (6) or oxygen (8)
    # and atom formal charge is zero.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8):
            return False, f"Extraneous atom found: {atom.GetSymbol()} (only C and O allowed)"
        if atom.GetFormalCharge() != 0:
            return False, f"Atom {atom.GetSymbol()} has non-zero formal charge"

    # 3. Count carbon atoms.
    atoms = mol.GetAtoms()
    c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, f"Too few carbon atoms ({c_count}); monosaccharides require 3–9 carbons"
    if c_count > 9:
        return False, f"Too many carbon atoms ({c_count}); parent monosaccharides are expected to have 3–9 carbons"
    
    # 4. Count hydroxyl groups using SMARTS "[OX2H]"
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(hydroxyl_smarts)
    num_oh = len(oh_matches)
    oh_ratio = num_oh / c_count

    # 5. Detect a cyclic sugar ring if present.
    ring_found = False
    ring_size = None
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            # Check if there is at least one oxygen in the ring.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                ring_found = True
                ring_size = len(ring)
                break

    # 6. Check for open-chain carbonyl signatures.
    # For terminal aldehyde: pattern "[H]C(=O)"
    aldehyde_smarts = Chem.MolFromSmarts("[H]C(=O)")
    has_terminal_aldehyde = mol.HasSubstructMatch(aldehyde_smarts)
    # For internal ketone: pattern "C(=O)C"
    ketone_smarts = Chem.MolFromSmarts("C(=O)C")
    has_internal_ketone = mol.HasSubstructMatch(ketone_smarts)

    # Decision branch:
    details = f"Carbons: {c_count}, Hydroxyls: {num_oh} (OH/C ratio: {oh_ratio:.2f})"
    if ring_found:
        # With a cyclic structure, expect at least 3 -OH groups and ratio >=0.5.
        if num_oh < 3:
            return False, f"Cyclic {ring_size}-membered sugar ring found but too few –OH groups ({num_oh}); expected at least 3. {details}"
        if oh_ratio < 0.5:
            return False, f"Cyclic sugar ring found but OH/C ratio ({oh_ratio:.2f}) is below 0.5. {details}"
        return True, f"Matches cyclic monosaccharide criteria: {details}. Contains a {ring_size}-membered ring with oxygen."
    else:
        # Open-chain: must contain a carbonyl group.
        if has_terminal_aldehyde:
            # For an aldose, most carbons (except the carbonyl) should be hydroxylated.
            expected_oh = c_count - 1
            if num_oh < expected_oh:
                return False, f"Open-chain aldose (terminal aldehyde detected) but only {num_oh} –OH groups (expected ~{expected_oh}). {details}"
            if oh_ratio < 0.5:
                return False, f"Open-chain aldose detected but OH/C ratio ({oh_ratio:.2f}) is too low. {details}"
            return True, f"Matches open-chain aldose criteria: {details}. Contains a terminal aldehyde."
        elif has_internal_ketone:
            # For a ketose, expect somewhat fewer OH groups (rough heuristic: expected ~C-2).
            expected_oh = c_count - 2
            if num_oh < expected_oh:
                return False, f"Open-chain ketose (internal ketone detected) but only {num_oh} –OH groups (expected ~{expected_oh}). {details}"
            if oh_ratio < 0.5:
                return False, f"Open-chain ketose detected but OH/C ratio ({oh_ratio:.2f}) is too low. {details}"
            return True, f"Matches open-chain ketose criteria: {details}. Contains an internal ketone."
        else:
            return False, "No cyclic sugar ring or recognizable open-chain carbonyl (aldehyde/ketone) detected."

# Example usage for testing
if __name__ == "__main__":
    test_examples = {
        # True positives (monosaccharides)
        "beta-ascarylopyranose": "C[C@@H]1O[C@H](O)[C@H](O)C[C@H]1O",
        "aldehydo-D-fucose": "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@@H](C)O",
        "alpha-D-idopyranose": "O1[C@@H]([C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)CO",
        "aldehydo-L-glucose": "[H]C(=O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)CO",
        "beta-D-gulose": "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O",
        "beta-L-sorbofuranose": "OC[C@@H]1O[C@@](O)(CO)[C@@H](O)[C@@H]1O",
        "D-fucopyranose": "C[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
        "L-altrofuranose": "O1[C@H]([C@H](O)[C@@H](O)C1O)[C@@H](O)CO",
        "alpha-L-galactose": "OC[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O",
        "2-dehydro-D-glucopyranose": "OC[C@H]1OC(O)C(=O)[C@@H](O)[C@@H]1O",
        "6-deoxy-beta-L-talopyranose": "O1[C@H]([C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O)C",
        "L-lyxofuranose": "O1[C@H]([C@@H](O)[C@@H](O)C1O)CO",
        "L-allofuranose": "[H][C@]1(OC(O)[C@@H](O)[C@H]1O)[C@@H](O)CO",
        "beta-D-fructuronic acid": "OC[C@@]1(O)O[C@@H]([C@@H](O)[C@@H]1O)C(O)=O",
        "2,3,4,5-tetrahydroxypentanal": "OCC(O)C(O)C(O)C=O",
        "D-galacto-hexodialdose": "O[C@@H](C=O)[C@@H](O)[C@@H](O)[C@H](O)C=O",
        "D-fructofuranuronic acid": "OCC1(O)O[C@@H]([C@@H](O)[C@@H]1O)C(O)=O",
        "aldehydo-L-xylose": "[H]C(=O)[C@@H](O)[C@H](O)[C@@H](O)CO",
        "aldehydo-L-iduronic acid": "O[C@@H](C=O)[C@@H](O)[C@H](O)[C@@H](O)C(O)=O",
        "D-altro-D-manno-Heptose": "O1C(C(O)C(O)C(O)C1O)C(O)CO",
        "alpha-D-talofuranose": "O1[C@H]([C@H](O)[C@H](O)[C@H]1O)[C@H](O)CO",
        "D-xylulose": "C(O)C(=O)[C@@H](O)[C@H](O)CO",
        "beta-D-glucofuranose": "O1[C@@]([C@H](O)[C@@H](O)[C@@H]1O)([C@H](O)CO)[H]",
        "alpha-D-mannose": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
        "beta-D-fructopyranose": "OC[C@@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",
        # False positives (should be rejected)
        "Diglycolic acid": "O(CC(O)=O)CC(O)=O",
        "alpha-(hydroxymethyl)serine zwitterion": "C(C([O-])=O)(CO)(CO)[NH3+]"
    }
    for name, smi in test_examples.items():
        result, msg = is_monosaccharide(smi)
        print(f"{name}: {result} | {msg}")