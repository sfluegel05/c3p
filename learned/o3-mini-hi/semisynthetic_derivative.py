"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: Semisynthetic derivative
Definition: Any organic molecular entity derived from a natural product by partial chemical synthesis.
This improved function uses a weighted scoring system rather than hard thresholds.
Heuristics considered:
  - Valid organic molecule (must contain carbon)
  - Molecular weight between 200 and 2000 Da (common for natural products)
  - Presence of a ring system (most natural products have rings)
  - Presence of chiral centers (often present in natural product scaffolds; not strictly required)
  - A moderate-to-high fraction of sp3 carbons (flat, fully aromatic molecules are penalized)
  - Excessive flexibility (too many rotatable bonds is penalized)
Scores are accumulated and if the overall score exceeds a threshold then the molecule is likely a semisynthetic derivative.
Note: These rules remain heuristic and can misclassify some compounds.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is likely to be a semisynthetic derivative based on its SMILES string.
    Instead of hard cutoffs, this revised function uses a weighted score.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if scoring suggests a semisynthetic derivative, False otherwise.
        str: Explanation of the scoring details and reasoning.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule is organic
    atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if 6 not in atom_nums:
        return False, "No carbon atoms found; unlikely to be organic"

    scoring_details = []
    score = 0

    # Check molecular weight: reward if within 200-2000 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if 200 <= mol_wt <= 2000:
        score += 1
        scoring_details.append(f"+1: Molecular weight {mol_wt:.1f} Da is within range")
    else:
        scoring_details.append(f"0: Molecular weight {mol_wt:.1f} Da is outside typical range")

    # Check for ring systems: most natural products have rings.
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings >= 1:
        score += 1
        scoring_details.append(f"+1: Found {num_rings} ring(s)")
    else:
        score -= 1
        scoring_details.append(" -1: No ring system found")

    # Check for chiral centers. Rather than strictly requiring chiral centers, we reward them if present.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) >= 1:
        score += 1
        scoring_details.append(f"+1: Found {len(chiral_centers)} chiral center(s)")
    else:
        # No chiral centers reduces confidence but is not an automatic rejection.
        score -= 0.5
        scoring_details.append(" -0.5: No chiral centers found")

    # Check fraction of sp3 carbons: low sp3 fraction may indicate a synthetic, flat structure.
    try:
        frac_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    except Exception:
        frac_sp3 = None
    if frac_sp3 is not None:
        if frac_sp3 >= 0.3:
            score += 1
            scoring_details.append(f"+1: Fraction of sp3 carbons is {frac_sp3:.2f}")
        else:
            score -= 0.5
            scoring_details.append(f" -0.5: Fraction of sp3 carbons is low ({frac_sp3:.2f})")
    else:
        scoring_details.append("0: Could not compute fraction of sp3 carbons")
    
    # Check number of rotatable bonds: too flexible molecules are less typical.
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds <= 30:
        score += 0.5
        scoring_details.append(f"+0.5: {rot_bonds} rotatable bond(s) within acceptable range")
    else:
        score -= 0.5
        scoring_details.append(f" -0.5: Too many rotatable bonds ({rot_bonds})")

    # Define a threshold for classification.
    threshold = 2.0  # Adjust threshold as needed for balance
    decision = score >= threshold

    # Build a detailed reasoning message.
    detail_msg = "; ".join(scoring_details) + f" => Total score: {score:.2f} (threshold: {threshold})"
    if decision:
        return True, ("Molecule is likely a semisynthetic derivative based on heuristics: " + detail_msg)
    else:
        return False, ("Molecule does not appear to be a semisynthetic derivative: " + detail_msg)

# Example usage:
# result, explanation = is_semisynthetic_derivative("CC[C@H]1OC(=O)[C@H](C)...")
# print(result, explanation)