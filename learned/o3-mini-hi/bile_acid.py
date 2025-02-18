"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: Hydroxy-5β-cholanic acids (bile acids)
Definition: Bile acids are members of the hydroxy-5β-cholanic acid family (or their glycine/taurine conjugates). 
They typically show a fused steroid nucleus built from one 5-membered ring and three 6-membered rings (roughly 17 carbon atoms total in the core), 
a carboxyl (or conjugated amide) group on a side chain (the acid carbon is not part of the core and is 3–7 bonds away), and at least three hydroxyl groups on the scaffold.
Note: This heuristic approach may need further refinement (for example, stereochemical tests to rule out 5α-configuration).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a (hydroxy-5β-cholanic) bile acid or its glycine/taurine conjugate.

    Heuristic tests:
      1. The SMILES string must parse properly.
      2. Contains exactly one acid-like moiety – either a free carboxylic acid or a glycine/taurine-conjugated acyl group
         (using SMARTS patterns, ensuring the acid group is not part of any ring).
      3. Must have at least three hydroxyl groups ([OX2H] groups).
      4. The molecular weight should be within the bile acid range (roughly 350–700 Da).
      5. Contains a fused steroid nucleus. Typically, this means one 5-membered ring and three 6-membered rings. 
         When only considering carbon atoms from these rings, a total count of roughly 16–18 is expected.
      6. The acid (or conjugate) group must be on a side chain – that is, its nearest distance (in number of bonds) to any nucleus carbon should be 3–7.
         
    Args:
       smiles (str): SMILES string for the molecule.
       
    Returns:
       (bool, str): Tuple with True and a positive reason if the molecule is classified as a bile acid,
                    or False and a reason for rejection.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Acid / Conjugate Group Checks ---
    # Free carboxylic acid group (not in a ring)
    acid_pattern = Chem.MolFromSmarts("[$([CX3](=O)[O;H,-]);!R]")
    # Glycine conjugate: acyl carbon bonded to -NCC(=O)[O;H,-]
    glycine_pattern = Chem.MolFromSmarts("[$([CX3](=O)NCC(=O)[O;H,-]);!R]")
    # Taurine conjugate: acyl carbon bonded to -NCCS(=O)(=O)[O;H,-]
    taurine_pattern = Chem.MolFromSmarts("[$([CX3](=O)NCCS(=O)(=O)[O;H,-]);!R]")
    
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    taurine_matches = mol.GetSubstructMatches(taurine_pattern)
    total_acid_matches = len(acid_matches) + len(glycine_matches) + len(taurine_matches)
    
    if total_acid_matches == 0:
        return False, "Missing a free acid or recognized conjugated acid group essential for cholanic acids"
    if total_acid_matches > 1:
        return False, f"Found {total_acid_matches} acid-like groups; expected exactly one."
    
    # --- Hydroxyl Group Count ---
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Too few hydroxyl substituents (found {len(hydroxyl_matches)}; expected at least 3)"
    
    # --- Molecular Weight Check ---
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (350 <= mol_wt <= 700):
        return False, f"Molecular weight {mol_wt:.1f} is out of the expected bile acid range (350–700 Da)"
    
    # --- Fused Steroid Nucleus Detection ---
    # Gather ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found; a steroid nucleus is expected"

    # Filter rings for common steroid sizes (5- and 6-membered rings)
    rings_5 = [ring for ring in rings if len(ring) == 5]
    rings_6 = [ring for ring in rings if len(ring) == 6]
    if len(rings_5) != 1 or len(rings_6) < 3:
        return False, f"Ring system does not meet steroid criteria (found {len(rings_5)} 5-membered and {len(rings_6)} 6-membered rings)"
    
    # Construct the set of atoms in our steroid nucleus (5- and 6-membered rings)
    steroid_core_atoms = set()
    for ring in rings:
        if len(ring) in (5, 6):
            steroid_core_atoms.update(ring)
    # Count only carbon atoms from the steroid core
    core_carbons = [atom for idx, atom in enumerate(mol.GetAtoms()) if idx in steroid_core_atoms and atom.GetAtomicNum() == 6]
    if not (16 <= len(core_carbons) <= 18):
        return False, f"Fused steroid core carbon count ({len(core_carbons)}) is not in the expected range (16–18)"
    
    # --- Side Chain Check for the Acid Group ---
    # Verify that the acid (or conjugate) group is on a side chain.
    # We take the acyl carbon (the first atom in any match) and verify that:
    #   (a) It is not part of the steroid nucleus.
    #   (b) Its shortest path (in bonds) to any nucleus carbon is between 3 and 7.
    acid_on_side_chain = False
    all_matches = list(acid_matches) + list(glycine_matches) + list(taurine_matches)
    for match in all_matches:
        acid_carbon_idx = match[0]
        if acid_carbon_idx in steroid_core_atoms:
            continue  # Acid carbon is part of the core, which is not acceptable.
        min_dist = float("inf")
        for core_idx in steroid_core_atoms:
            sp = Chem.GetShortestPath(mol, acid_carbon_idx, core_idx)
            if sp:
                dist = len(sp) - 1  # number of bonds in the path
                if dist < min_dist:
                    min_dist = dist
        if 3 <= min_dist <= 7:
            acid_on_side_chain = True
            break
    if not acid_on_side_chain:
        return False, "Acid/conjugate group is not situated on a side chain (expected separation of 3–7 bonds from the steroid core)"
    
    # If all tests pass, classify the molecule as a bile acid.
    return True, "Molecule displays a fused steroid nucleus (1 five-membered and at least 3 six-membered rings with ~17 carbons), " \
                  "a single, properly placed acid/conjugate group, and sufficient hydroxyl substituents typical for bile acids"

# --- Optional Testing Block ---
if __name__ == "__main__":
    # Testing with sample SMILES (true positives from provided examples)
    examples = [
        # (2S,3S,6R)-3-Hydroxy-2-methyl-6-((...)-heptanoic acid:
        "O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](C[C@H]3O)(C[C@H](O)CC4)[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CC[C@H](O)[C@H](C)C(O)=O)C)[H])[H])C",
        # 3alpha,11beta,12alpha-Trihydroxy-5beta-cholan-24-oic Acid:
        "O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(C[C@H](O)CC4)[H])C)([C@H]1O)[H])[H])(CC[C@@]2([C@@H](CCC(O)=O)C)[H])[H])C"
    ]
    for smi in examples:
        res, reason = is_bile_acid(smi)
        print(res, reason)