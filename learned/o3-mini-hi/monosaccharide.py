"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: Parent monosaccharide.
A parent monosaccharide is defined (heuristically) as either:
  (a) a polyhydroxy aldehyde (or ketone) in open‐chain form – having an open-chain aldehyde (or internal ketone)
      plus at least a majority of the carbons decorated with –OH groups (ratio ≥ 0.5) or
  (b) a cyclic sugar (5‐ or 6‐membered ring containing oxygen) with several –OH groups.
Only compounds with 3–9 carbon atoms and no extraneous heteroatoms (like P or S) are allowed.
Note: This is a heuristic; if the task is too difficult we may return (None, None).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a parent monosaccharide based on its SMILES string.
    Heuristics used:
      1. Must be a single fragment (no glycosidic bonds to other units).
      2. Carbon atom count must be between 3 and 9.
      3. Reject if containing extraneous heteroatoms such as phosphorus (P) or sulfur (S).
      4. Look for a sugar-like cyclic structure: a 5- or 6-membered ring that contains at least one oxygen.
         - If found, require at least 3 hydroxyl groups (SMARTS "[OX2H]") and a hydroxyl/carbon ratio ≥ 0.5.
      5. If no typical sugar ring is found, look for an open-chain carbonyl:
         - For aldoses, require a terminal aldehyde pattern ("[H]C(=O)").
         - For ketoses, require the presence of a carbonyl (pattern "[CX3]=[OX1]") but also a hydroxyl/carbon ratio ≥ 0.5.
      6. Reject molecules that have evidence of additional fragments (e.g. multiple disconnected parts) or aromatic carbonyls.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a parent monosaccharide, False otherwise.
        str: Reason for the classification decision.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return False, "Invalid SMILES string (exception during parsing)"
    if mol is None:
        return False, "Invalid SMILES string (could not parse molecule)"
    
    # 1. Check connectivity: must be a single fragment.
    frags = Chem.GetMolFrags(mol, asMols=False)
    if len(frags) > 1:
        return False, "Multiple fragments detected; likely not a single monosaccharide unit"
    
    # 2. Reject molecules containing extraneous heteroatoms: phosphorus (15) or sulfur (16)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in (15, 16):
            return False, f"Contains extraneous heteroatom {atom.GetSymbol()}; likely not a parent monosaccharide"
    
    # 3. Count carbons.
    atoms = mol.GetAtoms()
    c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, f"Too few carbon atoms ({c_count}); a monosaccharide must have at least 3 carbons"
    if c_count > 9:
        return False, f"Too many carbon atoms ({c_count}); parent monosaccharides are expected to have 3–9 carbons"
    
    # 4. Count hydroxyl groups.
    # Use a SMARTS for hydroxyl: oxygen with hydrogen attached.
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(hydroxyl_smarts)
    num_oh = len(oh_matches)
    
    # Calculate ratio of hydroxyl groups per carbon.
    oh_ratio = num_oh / c_count if c_count > 0 else 0
    
    # 5. Check for cyclic sugar rings:
    ring_found = False
    ring_type = None
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) in (5, 6):
            # Verify that at least one atom in the ring is oxygen.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                ring_found = True
                ring_type = f"{len(ring)}-membered"
                break
    
    # 6. Check for carbonyl group(s). We look for any [CX3]=[OX1].
    carbonyl_smarts = Chem.MolFromSmarts("[CX3]=[OX1]")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_smarts)
    
    # Also, for aldoses, we want a terminal aldehyde: pattern "[H]C(=O)"
    aldehyde_smarts = Chem.MolFromSmarts("[H]C(=O)")
    has_terminal_aldehyde = mol.HasSubstructMatch(aldehyde_smarts)
    
    # 7. Reject if a conjugated aromatic carbonyl is detected.
    aromatic_carbonyl_smarts = Chem.MolFromSmarts("[$([c][CX3]=[OX1])]")
    if mol.HasSubstructMatch(aromatic_carbonyl_smarts):
        return False, "Contains aromatic carbonyl; likely not a parent monosaccharide"
    
    # Decision making:
    reason_details = f"Carbons: {c_count}, Hydroxyls: {num_oh} (ratio {oh_ratio:.2f})"
    if ring_found:
        # With a typical sugar ring, even if the OH count is a bit low we may accept.
        if num_oh < 3:
            return False, f"Ring ({ring_type}) found but too few hydroxyl groups detected ({num_oh}); expect at least 3 for a cyclic sugar. {reason_details}"
        if oh_ratio < 0.5:
            return False, f"Ring found but the hydroxyl/carbon ratio ({oh_ratio:.2f}) is too low; expected ≥ 0.5. {reason_details}"
        extra = f"Contains a typical cyclic sugar ring ({ring_type})."
        return True, f"Matches monosaccharide criteria: {reason_details}. {extra}"
    else:
        # Open-chain case. We then require a proper carbonyl group.
        if has_terminal_aldehyde:
            # Likely an aldose.
            expected_oh = c_count - 1  # most aldoses have one carbonyl and remaining OH-attached carbons.
            if num_oh < expected_oh:
                return False, f"Open-chain aldose detected (terminal aldehyde) but insufficient hydroxyl groups ({num_oh} vs expected at least {expected_oh}). {reason_details}"
            if oh_ratio < 0.5:
                return False, f"Open-chain aldose detected but hydroxyl/carbon ratio ({oh_ratio:.2f}) is too low; expected at least 0.5. {reason_details}"
            extra = "Contains a terminal aldehyde indicating an open-chain aldose."
            return True, f"Matches monosaccharide criteria: {reason_details}. {extra}"
        elif has_carbonyl:
            # Possibly a ketose. In that case we require sufficient OH groups.
            expected_oh = c_count - 2  # ketoses typically have one less OH than aldoses.
            if num_oh < expected_oh:
                return False, f"Open-chain ketose detected (internal carbonyl) but insufficient hydroxyl groups ({num_oh} vs expected at least {expected_oh}). {reason_details}"
            if oh_ratio < 0.5:
                return False, f"Open-chain ketose detected but hydroxyl/carbon ratio ({oh_ratio:.2f}) is too low; expected at least 0.5. {reason_details}"
            extra = "Contains an internal carbonyl indicating an open-chain ketose."
            return True, f"Matches monosaccharide criteria: {reason_details}. {extra}"
        else:
            return False, "No cyclic sugar ring or open-chain carbonyl group detected; not a parent monosaccharide"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_examples = {
        # True positives:
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
        # False positives (should be rejected):
        "Diglycolic acid": "O(CC(O)=O)CC(O)=O",
        "alpha-(hydroxymethyl)serine zwitterion": "C(C([O-])=O)(CO)(CO)[NH3+]"
    }
    for name, smi in test_examples.items():
        result, msg = is_monosaccharide(smi)
        print(f"{name}: {result} | {msg}")