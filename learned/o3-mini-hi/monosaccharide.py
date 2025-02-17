"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: Parent monosaccharide.
A parent monosaccharide is defined heuristically as either:
  (a) A cyclic sugar unit – a 5- or 6-membered ring that contains exactly one oxygen
      and at most one exocyclic carbon (as in a CH2OH branch) with a high –OH density, or
  (b) An open‐chain polyhydroxy aldehyde or ketone (but not a carboxylic acid derivative)
      with a terminal aldehyde (or internal ketone) and –OH groups approaching (C-1) (or C-2).
Only molecules with only C and O atoms, zero formal charges, a single fragment,
and 3–9 carbons are acceptable.
Heuristics use an allowed tolerance (70% of the full expected number of –OH groups)
to catch deoxy sugars.
If the rule is too difficult to apply, the function may return (None, None).
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a parent monosaccharide based on its SMILES string.
    
    Heuristics applied:
      1. Must be a single-fragment molecule with only carbon and oxygen atoms (H implicit),
         and all atoms having formal charge zero.
      2. The number of carbon atoms must be between 3 and 9.
      3. Look for a cyclic sugar ring – a 5- or 6-membered ring that contains exactly one oxygen.
         If found, then the total number of carbons must be either equal to the number in the ring
         or one greater (to allow a CH2OH substituent). In this case the expected hydroxyl (–OH)
         count is (total carbons – 1) and the molecule is accepted if the actual count is at least 70%.
      4. If no proper cyclic ring is found, then look for open-chain carbonyl groups:
            a. Terminal aldehyde: SMARTS "[H]C(=O)[C]" – expected –OH count = (C-1)
            b. Internal ketone: SMARTS "C(=O)C" – expected –OH count = (C-2)
         Accept if the corresponding pattern exists and the –OH count is at least 70% of the
         expected value.
      5. Otherwise, the molecule is rejected.
      
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        (bool, str): True with a positive reason if classified as a parent monosaccharide,
                     or False plus a message indicating why it was rejected.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "SMILES parsing failed"
    
    # 1. Check for single fragment.
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        return False, "Multiple fragments detected; not a single monosaccharide unit"
    
    # 2. Check that all atoms are either carbon (6) or oxygen (8) with zero formal charge.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8):
            return False, f"Extraneous atom found: {atom.GetSymbol()} (only C and O allowed)"
        if atom.GetFormalCharge() != 0:
            return False, f"Atom {atom.GetSymbol()} has non-zero formal charge"
    
    # 3. Count carbon atoms.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    if c_count < 3:
        return False, f"Too few carbon atoms ({c_count}); monosaccharides require 3–9 carbons"
    if c_count > 9:
        return False, f"Too many carbon atoms ({c_count}); parent monosaccharides are expected to have 3–9 carbons"
    
    # 4. Count hydroxyl groups using the SMARTS pattern "[OX2H]".
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(hydroxyl_smarts)
    num_oh = len(oh_matches)
    
    details = f"Carbons: {c_count}, Hydroxyls: {num_oh}"
    
    # Function to test if observed hydroxyl count meets expected threshold.
    def meets_oh_threshold(expected):
        # Allow tolerance: require at least 70% of the expected –OH groups.
        return num_oh >= max(1, int(round(expected * 0.7)))
    
    # 5. Look for cyclic sugar rings.
    ring_info = mol.GetRingInfo().AtomRings()
    for ring in ring_info:
        ring_size = len(ring)
        # Check if ring is 5- or 6-membered.
        if ring_size in (5, 6):
            # Check that exactly one oxygen is in the ring.
            oxygen_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxygen_in_ring != 1:
                continue
            # Determine how many carbons are in the ring.
            ring_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            # For a "pure" cyclic monosaccharide, total carbons should be equal 
            # to the ring carbons (if no exocyclic substituent) or one greater (allowing a CH2OH group).
            if c_count not in (ring_carbons, ring_carbons + 1):
                continue
            # Expected –OH groups: for a full sugar, typically (total carbons - 1).
            expected_oh = c_count - 1
            if not meets_oh_threshold(expected_oh):
                return False, f"Cyclic {ring_size}-membered ring found but too few –OH groups ({num_oh}); expected ~{expected_oh}. {details}"
            return True, f"Matches cyclic monosaccharide criteria: {details} (OH/C ratio: {num_oh/c_count:.2f}). Contains a {ring_size}-membered ring with one oxygen."
    
    # 6. No acceptable cyclic ring found so look for open-chain patterns.
    # Define SMARTS for terminal aldehyde and internal ketone.
    aldehyde_smarts = Chem.MolFromSmarts("[H]C(=O)[C]")
    ketone_smarts = Chem.MolFromSmarts("C(=O)C")
    has_ald = mol.HasSubstructMatch(aldehyde_smarts)
    has_ket = mol.HasSubstructMatch(ketone_smarts)
    
    if has_ald:
        # For an aldose, expect roughly (C-1) hydroxyl groups.
        expected_oh = c_count - 1
        if not meets_oh_threshold(expected_oh):
            return False, f"Open-chain aldose (terminal aldehyde detected) but insufficient –OH groups ({num_oh} vs expected ~{expected_oh}). {details}"
        return True, f"Matches open-chain aldose criteria: {details} (OH/C ratio: {num_oh/c_count:.2f}). Contains terminal aldehyde."
    elif has_ket:
        # For a ketose, expect roughly (C-2) –OH groups.
        expected_oh = c_count - 2
        if not meets_oh_threshold(expected_oh):
            return False, f"Open-chain ketose (internal ketone detected) but insufficient –OH groups ({num_oh} vs expected ~{expected_oh}). {details}"
        return True, f"Matches open-chain ketose criteria: {details} (OH/C ratio: {num_oh/c_count:.2f}). Contains internal ketone."
    else:
        return False, "No cyclic sugar ring or recognizable open-chain carbonyl (aldehyde/ketone) detected."

# Example usage for testing.
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