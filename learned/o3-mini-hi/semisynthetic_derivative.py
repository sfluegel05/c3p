"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: Semisynthetic derivative
Definition: Any organic molecular entity derived from a natural product by partial chemical synthesis.

This version uses a weighted scoring system with subtle adjustments:
  • Molecular weight between ~200 and 2000 Da gives a +1 bonus. However, for molecules with MW < 250 Da we now subtract an extra 0.5.
  • At least one ring is required (+1) or else −1.
  • Chiral centers contribute:
        – For molecules with MW <300 Da: if at least one chiral center is found, add +0.5.
        – Otherwise (MW ≥300 Da): if there are ≥2 chiral centers add +1, else subtract 0.5.
  • A high fraction of sp3 carbons (≥0.3) gives +1; else −0.5.
  • Rotatable bonds:
        – ≤8 bonds: +1,
        – ≤15 bonds: +0.5,
        – >15: −1.0.
  • Heteroaromatic bonus: if any aromatic ring (with at least 5 atoms) contains at least one non‐carbon, add +0.5.
  • Chiral density (ratio of chiral centers to heavy atoms) is rewarded only if ≥0.15, giving +1.
The overall score must meet a tougher threshold (set here to 4.0) to be classified as a semisynthetic derivative.
This tuning was motivated by our past outcomes (where many false positives scored around 2.5–4.0) and aims to reduce mis‐classifications.
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
        str: A detailed explanation of the scoring.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic sanity: molecule should be organic (contain carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atoms found; unlikely to be organic"

    scoring_details = []  # to record the reasoning steps
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
        score -= 0.5  # increased penalty compared to previous version
        scoring_details.append(" -0.5: Molecular weight below 250 Da")

    # 2. Ring systems: most natural products have at least one ring.
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings >= 1:
        score += 1.0
        scoring_details.append(f"+1: Found {num_rings} ring(s)")
    else:
        score -= 1.0
        scoring_details.append(" -1: No ring system found")
    
    # 3. Chiral centers: use RDKit's chiral center finder.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral = len(chiral_centers)
    # For small molecules (MW <300), do not penalize lacking chiral centers.
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
    
    # 4. Fraction of sp3 carbons. A high fraction suggests a natural product–like 3D scaffold.
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
    
    # 5. Rotatable bonds: too much flexibility reduces the likelihood.
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds <= 8:
        rot_score = 1.0
        scoring_details.append(f"+1: {rot_bonds} rotatable bond(s) (very rigid)")
    elif rot_bonds <= 15:
        rot_score = 0.5
        scoring_details.append(f"+0.5: {rot_bonds} rotatable bond(s) within acceptable range")
    else:
        rot_score = -1.0  # harsher penalty for very flexible molecules
        scoring_details.append(f" -1.0: Too many rotatable bonds ({rot_bonds})")
    score += rot_score

    # 6. Heteroaromatic rings: if any aromatic ring (with at least 5 atoms) contains a non–carbon, add bonus.
    bonus_given = False
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) < 5:
            continue  # ignore very small rings
        # Check that every atom in the ring is aromatic and that at least one is non–carbon.
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        if all(atom.GetIsAromatic() for atom in atoms_in_ring):
            if any(atom.GetAtomicNum() != 6 for atom in atoms_in_ring):
                score += 0.5
                scoring_details.append("+0.5: Found heteroaromatic ring")
                bonus_given = True
                break
    if not bonus_given:
        scoring_details.append("0: No heteroaromatic ring bonus")

    # 7. Chiral density bonus: ratio of chiral centers to heavy atoms.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if heavy_atoms:
        chiral_density = num_chiral / len(heavy_atoms)
        # Reward only if chiral density is at least 0.15.
        if chiral_density >= 0.15:
            score += 1.0
            scoring_details.append(f"+1: Chiral density is {chiral_density:.2f}")
        else:
            scoring_details.append(f"0: Chiral density is low ({chiral_density:.2f})")
    else:
        scoring_details.append("0: No heavy atoms to compute chiral density")
    
    # Set a stricter classification threshold.
    threshold = 4.0
    decision = score >= threshold

    # Compile the detailed reasoning.
    detail_msg = "; ".join(scoring_details) + f" => Total score: {score:.2f} (threshold: {threshold})"
    if decision:
        return True, "Molecule is likely a semisynthetic derivative based on heuristics: " + detail_msg
    else:
        return False, "Molecule does not appear to be a semisynthetic derivative: " + detail_msg

# Example usage (testing on a few molecules; in practice, the test set is extensive)
if __name__ == '__main__':
    test_smiles = {
        "roxithromycin": "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@@H](C)C(=NOCOCCOC)[C@H](C)[C@@H](O)[C@]1(C)O",
        "TX-1875": "C(\\[C@H]([C@]1([C@H](/C=C/CCC[C@H](C)O)C[C@@H](C1)O)[H])O)=C/C(=O)O",
        "N-ethylharmine": "CCn1c2cc(OC)ccc2c2ccnc(C)c12",
        "Conoideochromane B": "O1C2=C(O)C=C(COC)C=C2C[C@@H](C1(C)C)O"
    }
    for name, smi in test_smiles.items():
        result, explanation = is_semisynthetic_derivative(smi)
        print(f"{name}: {result}\n  Explanation: {explanation}\n")