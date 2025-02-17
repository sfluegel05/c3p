"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: Semisynthetic derivative
Definition: Any organic molecular entity derived from a natural product by partial chemical synthesis.

This function uses a weighted scoring system. In our reasoning we note that semisynthetic derivatives
often have:
  • a molecular weight between about 200 and 2000 Da (with a bonus if above ~250 Da),
  • at least one ring system,
  • many natural‐product scaffolds have chiral centers—but some semisynthetic alkaloids (eg N‑ethylharmine)
    may lack any; therefore if MW<300 we do not penalize missing chiral centers,
  • a moderate-to-high fraction of sp3 carbons (>0.3),
  • not too “floppy” a structure (we reward molecules with few rotatable bonds),
  • and sometimes the presence of a heteroaromatic ring (eg when a nitrogen is present in an aromatic ring)
    is a hint.
Because these rules are heuristic, no single cutoff is perfect.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines whether the given SMILES is consistent with a semisynthetic derivative.
    Uses a weighted scoring scheme with several loosely tuned criteria.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the overall score equals or exceeds the threshold.
        str: Explanation of the scoring (or, if SMILES cannot be parsed, an error message).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic sanity: molecule should be organic (contain carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atoms found; unlikely to be organic"

    scoring_details = []
    score = 0.0

    # 1. Molecular Weight: reward if between 200 and 2000 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if 200 <= mol_wt <= 2000:
        score += 1.0
        scoring_details.append(f"+1: Molecular weight {mol_wt:.1f} Da is within range")
    else:
        scoring_details.append(f"0: Molecular weight {mol_wt:.1f} Da is outside typical range")
    # Extra penalty for very low weight (<250 Da)
    if mol_wt < 250:
        score -= 0.5
        scoring_details.append(" -0.5: Molecular weight below 250 Da")

    # 2. Ring systems: most natural products have rings.
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings >= 1:
        score += 1.0
        scoring_details.append(f"+1: Found {num_rings} ring(s)")
    else:
        score -= 1.0
        scoring_details.append(" -1: No ring system found")
    
    # 3. Chiral centers:
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral = len(chiral_centers)
    # If the molecule is small (<300 Da) we do not harshly penalize missing chiral centers.
    if mol_wt < 300:
        if num_chiral >= 1:
            score += 0.5
            scoring_details.append(f"+0.5: Found {num_chiral} chiral center(s) (small molecule)")
        else:
            scoring_details.append("0: No chiral centers found (small molecule; not penalized)")
    else:
        if num_chiral >= 2:
            score += 1.0
            scoring_details.append(f"+1: Found {num_chiral} chiral center(s)")
        else:
            score -= 0.5
            scoring_details.append(f" -0.5: Low number of chiral centers ({num_chiral}) for a larger molecule")
    
    # 4. Fraction of sp3 carbons: high fraction (>=0.3) suggests a natural product–like 3D scaffold.
    try:
        frac_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    except Exception:
        frac_sp3 = None
    if frac_sp3 is not None:
        if frac_sp3 >= 0.3:
            score += 1.0
            scoring_details.append(f"+1: Fraction of sp3 carbons is {frac_sp3:.2f}")
        else:
            score -= 0.5
            scoring_details.append(f" -0.5: Fraction of sp3 carbons is low ({frac_sp3:.2f})")
    else:
        scoring_details.append("0: Could not compute fraction of sp3 carbons")
    
    # 5. Rotatable bonds: too much flexibility is less typical.
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds <= 8:
        rot_score = 1.0
        scoring_details.append(f"+1: {rot_bonds} rotatable bond(s) (very rigid)")
    elif rot_bonds <= 15:
        rot_score = 0.5
        scoring_details.append(f"+0.5: {rot_bonds} rotatable bond(s) within acceptable range")
    else:
        rot_score = -0.5
        scoring_details.append(f" -0.5: Too many rotatable bonds ({rot_bonds})")
    score += rot_score

    # 6. Bonus from heteroaromatic rings. Many semisynthetic alkaloids show an aromatic ring
    # with a non‐carbon atom (eg nitrogen). (This rule is applied only if at least one ring contains
    # a non‐carbon atom.)
    ring_info = mol.GetRingInfo()
    bonus_given = False
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring):
            score += 0.5
            scoring_details.append("+0.5: Found heteroaromatic ring")
            bonus_given = True
            break
    if not bonus_given:
        scoring_details.append("0: No heteroaromatic ring bonus")
        
    # (Optional) 7. Bonus from chiral density: ratio of chiral centers to heavy atoms.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if heavy_atoms:
        chiral_density = num_chiral / len(heavy_atoms)
        if chiral_density >= 0.1:
            score += 0.5
            scoring_details.append(f"+0.5: Chiral density is {chiral_density:.2f}")
        else:
            scoring_details.append(f"0: Chiral density is low ({chiral_density:.2f})")
    
    # Set a classification threshold.
    # (This threshold is “tunable” – here we use 2.5 as a balance point.)
    threshold = 2.5
    decision = score >= threshold

    # Build the detailed reasoning string.
    detail_msg = "; ".join(scoring_details) + f" => Total score: {score:.2f} (threshold: {threshold})"
    if decision:
        return True, "Molecule is likely a semisynthetic derivative based on heuristics: " + detail_msg
    else:
        return False, "Molecule does not appear to be a semisynthetic derivative: " + detail_msg

# Example usage:
if __name__ == '__main__':
    # Test a few of the example SMILES (only a few are shown; the full test set is extensive)
    test_smiles = {
        "roxithromycin": "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@@H](C)C(=NOCOCCOC)[C@H](C)[C@@H](O)[C@]1(C)O",
        "N-ethylharmine": "CCn1c2cc(OC)ccc2c2ccnc(C)c12",
        "Conoideochromane B": "O1C2=C(O)C=C(COC)C=C2C[C@@H](C1(C)C)O"
    }
    for name, smi in test_smiles.items():
        result, explanation = is_semisynthetic_derivative(smi)
        print(f"{name}: {result}\n  Explanation: {explanation}\n")