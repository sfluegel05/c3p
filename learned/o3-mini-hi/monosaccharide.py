"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: Parent monosaccharide

A parent monosaccharide is defined (heuristically) as either:
  (a) A cyclic sugar unit – a 5- or 6-membered ring that contains exactly one oxygen
      and at most one exocyclic carbon (to allow a CH2OH branch) with a high –OH density.
  (b) An open‐chain polyhydroxy aldehyde (or ketone) having a terminal aldehyde 
      (or an internal ketone) and –OH groups numbering near (total carbons –1) or (total carbons –2).
Additional requirements:
  • Only carbon and oxygen atoms are allowed (hydrogens implicit), and all atoms must have zero formal charge.
  • The molecule must be a single fragment.
  • The number of carbons must be between 3 and 9.
  • To catch deoxy sugars, a tolerance is allowed:
      – For cyclic sugars we require at least 60% of the ideal –OH number (ideal = total carbons – 1).
      – For open‐chain sugars we require at least 70% of the expected –OH number.
If the rules are hard to apply or ambiguous, the function may return (None, None).
"""

from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a parent monosaccharide based on its SMILES string.
    
    Args:
      smiles (str): SMILES string of the molecule
      
    Returns:
      (bool, str): (True, reason) if the molecule is classified as a parent monosaccharide,
                   (False, reason) if rejected, or (None, None) if the rule is too hard to apply.
    """
    # Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "SMILES parsing failed"
    
    # 1. Check that the molecule is a single fragment.
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        return False, "Multiple fragments detected; not a single monosaccharide unit"
        
    # 2. Ensure that every atom is either carbon or oxygen and that all atoms have formal charge zero.
    for atom in mol.GetAtoms():
        atomic = atom.GetAtomicNum()
        if atomic not in (6, 8):
            return False, f"Extraneous atom found: {atom.GetSymbol()} (only C and O allowed)"
        if atom.GetFormalCharge() != 0:
            return False, f"Atom {atom.GetSymbol()} has non-zero formal charge"
            
    # 3. Count carbons.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbon_atoms)
    if c_count < 3:
        return False, f"Too few carbon atoms ({c_count}); monosaccharides require 3–9 carbons"
    if c_count > 9:
        return False, f"Too many carbon atoms ({c_count}); parent monosaccharides are expected to have 3–9 carbons"
    
    # 4. Count hydroxyl groups using the SMARTS "[OX2H]".
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    num_oh = len(oh_matches)
    
    details = f"Carbons: {c_count}, Hydroxyls: {num_oh}, OH/C ratio: {num_oh/c_count:.2f}"
    
    # Helper function to decide if the number of –OH groups meets expected threshold.
    def meets_threshold(observed, expected, tol):
        # tol is the fraction required; use max(1,rounded_value)
        required = max(1, int(round(expected * tol)))
        return observed >= required
    
    # 5. Check for a cyclic sugar ring.
    ring_info = mol.GetRingInfo().AtomRings()
    for ring in ring_info:
        ring_size = len(ring)
        if ring_size not in (5, 6):
            continue  # only consider rings of size 5 or 6
        # Count oxygens in the ring.
        oxygens_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxygens_in_ring != 1:
            continue  # must have exactly one oxygen in the ring
        # Count carbons in the ring.
        ring_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        # Allow either a “pure” ring (c_count == ring_carbons)
        # or a ring with one extra carbon (e.g. a CH2OH substituent: c_count == ring_carbons + 1)
        if c_count not in (ring_carbons, ring_carbons + 1):
            continue
        # For a cyclic sugar, an ideal monosaccharide should have (total carbons – 1) –OH groups.
        expected_oh = c_count - 1
        # Use a 60% threshold here (to allow for deoxy variants).
        if meets_threshold(num_oh, expected_oh, 0.60):
            return True, f"Matches cyclic monosaccharide criteria: {details}. Contains a {ring_size}-membered ring with one oxygen."
        else:
            return False, (f"Cyclic {ring_size}-membered ring found but too few –OH groups "
                           f"({num_oh} vs expected ~{expected_oh} at 60% threshold). {details}")
    
    # 6. If no acceptable cyclic ring is found, try open-chain criteria.
    # Look for a terminal aldehyde or internal ketone.
    aldehyde_pattern = Chem.MolFromSmarts("[H]C(=O)[C]")
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    is_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    is_ketone = mol.HasSubstructMatch(ketone_pattern)
    
    if is_aldehyde:
        # For an aldose, expect about c_count - 1 –OH groups.
        expected_oh = c_count - 1
        if meets_threshold(num_oh, expected_oh, 0.70):
            return True, (f"Matches open-chain aldose criteria: {details} "
                          f"(expected –OH ~{expected_oh} at 70% threshold). Contains terminal aldehyde.")
        else:
            return False, (f"Open-chain aldose (terminal aldehyde detected) but insufficient –OH groups "
                           f"({num_oh} vs expected ~{expected_oh} at 70% threshold). {details}")
    elif is_ketone:
        # For a ketose, expect roughly (c_count - 2) –OH groups.
        expected_oh = c_count - 2
        if meets_threshold(num_oh, expected_oh, 0.70):
            return True, (f"Matches open-chain ketose criteria: {details} "
                          f"(expected –OH ~{expected_oh} at 70% threshold). Contains internal ketone.")
        else:
            return False, (f"Open-chain ketose (internal ketone detected) but insufficient –OH groups "
                           f"({num_oh} vs expected ~{expected_oh} at 70% threshold). {details}")
    else:
        return False, "No cyclic sugar ring or recognizable open-chain carbonyl (aldehyde/ketone) detected."

# Example usage (for testing):
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