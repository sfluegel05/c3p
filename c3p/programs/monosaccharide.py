"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: Monosaccharide (parent monosaccharides are defined as polyhydroxy 
aldehydes/ketones with three to about nine carbons and not covalently linked to other sugar units)
Based on a set of heuristics improved to pick up deoxy sugars and reject substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a parent monosaccharide based on its SMILES string.
    We use several heuristics:
      1. The molecule must be a single fragment (no glycosidic connection).
      2. The carbon count should lie between 3 and 9.
      3. The number of oxygens should be in a plausible range relative to the carbons
         (allowed range: from (C - 2) to (C + 2)).
      4. The molecule should have several hydroxyl groups (using the SMARTS "[OX2H]").
      5. The molecule should either have (a) a sugar ring – a 5- or 6-membered ring that contains at least one oxygen
         (typical for pyranoses or furanoses) or (b) an open-chain carbonyl (aldehyde or ketone).
      6. The molecule should not contain extraneous heteroatoms such as phosphorus (P) or sulfur (S)
         which are common in phosphorylated/sulfonated sugars or glycosides.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a parent monosaccharide.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check connectivity: must be a single fragment.
    frags = Chem.GetMolFrags(mol, asMols=False)
    if len(frags) > 1:
        return False, "Multiple fragments detected; likely not a single monosaccharide unit"
    
    # Reject if extraneous heteroatoms (P or S) are present.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            return False, "Contains phosphorus; likely not a parent monosaccharide"
        if atom.GetAtomicNum() == 16:
            return False, "Contains sulfur; likely not a parent monosaccharide"
    
    # Count carbon and oxygen atoms.
    atoms = mol.GetAtoms()
    c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
    
    if c_count < 3:
        return False, f"Too few carbon atoms ({c_count}); a monosaccharide must have at least 3 carbons"
    if c_count > 9:
        return False, f"Too many carbon atoms ({c_count}); parent monosaccharides are expected to have 3–9 carbons"
    
    # Allow oxygen counts roughly within [c_count - 2, c_count + 2]
    lower_bound = max(c_count - 2, 1)
    upper_bound = c_count + 2
    if not (lower_bound <= o_count <= upper_bound):
        return False, f"Oxygen count ({o_count}) is not within the expected range ({lower_bound}–{upper_bound}) for {c_count} carbons"
    
    # Count hydroxyl groups using SMARTS "[OX2H]":
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(hydroxyl_smarts)
    if len(oh_matches) < max(c_count - 2, 2):
        return False, f"Too few hydroxyl groups detected ({len(oh_matches)}); expect several -OH groups on a monosaccharide"
    
    # Check for ring structure typical for sugars.
    ring_found = False
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) in (5, 6):  # typical furanose or pyranose ring
            # Check if at least one atom in ring is oxygen.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                ring_found = True
                break

    # Check for a carbonyl group (open-chain aldehyde/ketone).
    carbonyl_smarts = Chem.MolFromSmarts("[CX3]=[OX1]")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_smarts)
    
    # For a monosaccharide, if it is not cyclic then a carbonyl should be present.
    if not ring_found and not has_carbonyl:
        return False, "Neither a typical sugar ring (5- or 6-membered with an oxygen) nor a carbonyl group found; not a parent monosaccharide"
    
    # Construct a reason string.
    extra = ""
    if ring_found:
        extra = "Contains a typical sugar ring."
    elif not has_carbonyl:
        extra = "Might be in cyclic form with a latent (hemiacetal) carbonyl."
    else:
        extra = "Contains an open-chain carbonyl."
    
    reason = (f"Matches monosaccharide criteria: {c_count} carbons, {o_count} oxygens, "
              f"{len(oh_matches)} hydroxyl groups. {extra}")
    return True, reason

# Example usage:
if __name__ == "__main__":
    # Try some examples:
    examples = {
        "aldehydo-L-glucose": "[H]C(=O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)CO",
        "beta-ascarylopyranose": "C[C@@H]1O[C@H](O)[C@H](O)C[C@H]1O",
        "D-galactosamine 6-phosphate": "N[C@H]1C(O)O[C@H](COP(O)(O)=O)[C@H](O)[C@@H]1O",
    }
    for name, smi in examples.items():
        result, msg = is_monosaccharide(smi)
        print(f"{name}: {result} | {msg}")