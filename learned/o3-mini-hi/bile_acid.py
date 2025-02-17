"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: Hydroxy-5β-cholanic acids (bile acids)
Definition: Bile acids are members of the hydroxy-5β-cholanic acid family (or their glycine/taurine conjugates). They typically show a fused steroid nucleus built from one 5-membered ring and three 6-membered rings (whose carbon count is roughly 17), a carboxyl (or conjugated amide) group on an aliphatic side‐chain (the acid carbon must be located away from the steroid core by at least 3 bonds but not too long – roughly 3–7 bonds), and at least three hydroxyl groups on the scaffold.
Note: This is only a heuristic approach and further refinement (for example, specific stereochemistry to rule out 5α‐configuration) may be needed.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a (hydroxy-5β-cholanic) bile acid or its glycine/taurine conjugate.
    
    Tests (heuristic):
      1. Molecule is parsed properly.
      2. Contains exactly one acid-like moiety – either a free carboxylic acid or a glycine/taurine-conjugated acyl group.
         (SMARTS patterns are used; note we require that the acid group is not in a ring.)
      3. Shows at least three hydroxyl groups ([OX2H]).
      4. Has a molecular weight roughly in the bile acid range (350–700 Da).
      5. Contains a fused steroid nucleus defined from rings that are exactly one 5-membered ring and three 6-membered rings.
         Moreover, if only the carbons in these rings are considered, the count should be roughly 16–18 (typically ~17 for cholanic acids).
      6. The acid (or conjugate) group is on a side chain – i.e. its acyl carbon is not part of the nucleus and its shortest path
         (in bonds) to any nucleus carbon is at least 3 but not excessively long (3–7).
         
    Args:
       smiles (str): SMILES string for the molecule.
       
    Returns:
       (bool, str): Tuple with True and a positive reason if the molecule is classified as a bile acid, 
                    or False and a reason for rejection.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Acid/Conjugate group pattern(s) ---
    # Free carboxylic acid – ensure it is not in a ring.
    acid_pattern = Chem.MolFromSmarts("[$([CX3](=O)[O;H,-]);!R]")
    # Glycine conjugate: amide where acyl carbon is bonded to –NCC(=O)[O;H,-]
    glycine_pattern = Chem.MolFromSmarts("[$([CX3](=O)NCC(=O)[O;H,-]);!R]")
    # Taurine conjugate: acyl carbon bonded to –NCCS(=O)(=O)[O;H,-]
    taurine_pattern = Chem.MolFromSmarts("[$([CX3](=O)NCCS(=O)(=O)[O;H,-]);!R]")
    
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    taurine_matches = mol.GetSubstructMatches(taurine_pattern)
    total_acid_matches = len(acid_matches) + len(glycine_matches) + len(taurine_matches)
    if total_acid_matches == 0:
        return False, "Missing free acid or recognized conjugated acid group essential for cholanic acids"
    if total_acid_matches > 1:
        # If more than one acid-like group detected, we consider this ambiguous
        return False, f"Found {total_acid_matches} acid-like groups; expected exactly one."
    
    # --- Hydroxyl count ---
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Too few hydroxyl substituents (found {len(hydroxyl_matches)}; expect at least 3)"
    
    # --- Molecular Weight ---
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (350 <= mol_wt <= 700):
        return False, f"Molecular weight {mol_wt:.1f} out of expected bile acid range (350-700 Da)"
    
    # --- Fused Steroid Nucleus ---
    # Get all rings in the molecule.
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized():
        return False, "No ring information available"
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found; expected a steroid core"
    
    # Filter rings that are 5- or 6-membered.
    rings_5 = [ring for ring in rings if len(ring) == 5]
    rings_6 = [ring for ring in rings if len(ring) == 6]
    if len(rings_5) != 1 or len(rings_6) < 3:
        return False, f"Ring system does not meet steroid criteria (found {len(rings_5)} 5-membered and {len(rings_6)} 6-membered rings)"
    
    # Assemble a fused steroid nucleus from those rings.
    steroid_core_atoms = set()
    for ring in rings:
        if len(ring) in (5, 6):
            steroid_core_atoms.update(ring)
    # Count only carbon atoms in the core.
    core_carbons = [atom for idx, atom in enumerate(mol.GetAtoms()) if idx in steroid_core_atoms and atom.GetAtomicNum() == 6]
    if not (16 <= len(core_carbons) <= 18):
        return False, f"Fused steroid core carbon count ({len(core_carbons)}) not in expected range (16–18)"
    
    # --- Side-chain check ---
    # The acid (or conjugate) group should be on a side chain.
    # For our purposes, we take the first atom of the acid SMARTS match (the acyl carbon)
    # and require that it is NOT in the steroid core AND that its shortest path distance to any core atom is at least 3
    # and not overly long (we assume 3–7 bonds).
    all_matches = list(acid_matches) + list(glycine_matches) + list(taurine_matches)
    acid_on_side_chain = False
    for match in all_matches:
        acid_carbon_idx = match[0]
        if acid_carbon_idx in steroid_core_atoms:
            continue  # acid carbon is part of the core, which is not acceptable
        min_dist = float("inf")
        for core_idx in steroid_core_atoms:
            sp = Chem.GetShortestPath(mol, acid_carbon_idx, core_idx)
            if sp:
                dist = len(sp) - 1  # number of bonds
                if dist < min_dist:
                    min_dist = dist
        if 3 <= min_dist <= 7:
            acid_on_side_chain = True
            break
    if not acid_on_side_chain:
        return False, "Acid/conjugate group appears not to be on a proper side chain (expected separation of 3–7 bonds from steroid core)"
    
    # If all tests pass then we judge the molecule as a bile acid.
    return True, "Molecule shows a fused steroid nucleus (1 five-membered and 3 six-membered rings with ~17 carbons) " \
                  "with a single, properly placed acid/conjugate group and at least three hydroxyls typical for bile acids"


# --- Optional testing block ---
if __name__ == "__main__":
    # Some example SMILES from the true positive list:
    examples = [
        # (2S,3S,6R)-3-Hydroxy-2-methyl-6-((...)-heptanoic acid:
        "O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](C[C@H]3O)(C[C@H](O)CC4)[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CC[C@H](O)[C@H](C)C(O)=O)C)[H])[H])C",
        # 3alpha,11beta,12alpha-Trihydroxy-5beta-cholan-24-oic Acid:
        "O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(C[C@H](O)CC4)[H])C)([C@H]1O)[H])[H])(CC[C@@]2([C@@H](CCC(O)=O)C)[H])[H])C"
    ]
    for smi in examples:
        res, reason = is_bile_acid(smi)
        print(res, reason)