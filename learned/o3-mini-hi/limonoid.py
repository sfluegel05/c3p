"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid class – tetranortriterpenoids that are highly oxygenated and either contain 
or are derived from a precursor with a 4,4,8-trimethyl-17-furanylsteroid skeleton.
Heuristic criteria (each can contribute to a total score):
  1. Molecular weight between 350 and 750 Da. (1.0 point)
  2. Polycyclic structure:
       - At least 4 rings yields 1.0 point.
       - 6 or more rings give a bonus of 0.5.
  3. Carbon count between 23 and 40. (1.0 point)
  4. Oxygenation: if the oxygen-to-carbon ratio is at least 0.17 but not too high (≤0.55), add 1.0 point.
  5. Furan ring presence (SMARTS "c1occc1"). (1.0 point)
Additionally, if a sugar moiety is detected (via a pyranose‐like SMARTS), 
we subtract 1.5 points because such carbohydrate portions are typical of glycosides.
A total (adjusted) score of ≥4.0 (if a furan is present) or ≥3.5 (if furan is absent)
is taken as evidence to classify a molecule as a limonoid.
Note:
  – Some limonoids (e.g. defurano compounds) lose the furan bonus.
  – This heuristic is imperfect but is tuned to improve our F1 performance by reducing false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines whether a molecule belongs to the limonoid class based on heuristic criteria.
    
    The criteria (with associated points) are:
      1. MW between 350 and 750 Da (1.0 point).
      2. Rings: at least 4 rings (1.0 point), with a +0.5 bonus if the rings count is 6 or more.
      3. Carbon count between 23 and 40 (1.0 point).
      4. Oxygenation: if the oxygen-to-carbon (O/C) ratio is between 0.17 and 0.55 (inclusive), add 1.0 point.
      5. Presence of a furan ring (SMARTS "c1occc1") adds 1.0 point.
      
    Additionally, if a sugar (pyranose-type) substructure is detected, subtract 1.5 points.
    
    The final classification uses adjusted thresholds:
      - If a furan ring is found, a total score of >= 4.0 is classified as a limonoid.
      - If NO furan ring is found, a lower threshold of >= 3.5 is used.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple of classification (True if limonoid, False otherwise) and a reason text.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    reasons = []
    total_score = 0.0  # initialize score

    # Criterion 1: Molecular weight between 350 and 750 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if 350 <= mol_wt <= 750:
        total_score += 1.0
        reasons.append(f"MW {mol_wt:.1f} Da is within expected range (350-750).")
    else:
        reasons.append(f"MW {mol_wt:.1f} Da is outside expected range (350-750).")
    
    # Criterion 2: Count rings (polycyclic structure).
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings >= 4:
        total_score += 1.0
        reasons.append(f"Number of rings ({num_rings}) is sufficient (>=4).")
        if num_rings >= 6:
            total_score += 0.5
            reasons.append("Bonus: Rings count is high (>=6), adding extra 0.5.")
    else:
        reasons.append(f"Only {num_rings} rings found (expected at least 4).")
    
    # Criterion 3: Carbon count between 23 and 40.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if 23 <= carbon_count <= 40:
        total_score += 1.0
        reasons.append(f"Carbon count ({carbon_count}) is within expected range (23-40).")
    else:
        reasons.append(f"Carbon count ({carbon_count}) is outside expected range (23-40).")
    
    # Criterion 4: Oxygenation level (Oxygen-to-Carbon ratio)
    # We now add a point only if 0.17 <= O/C ratio <= 0.55.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    o_c_ratio = oxygen_count / carbon_count if carbon_count > 0 else 0
    if 0.17 <= o_c_ratio <= 0.55:
        total_score += 1.0
        reasons.append(f"O/C ratio ({o_c_ratio:.2f}) is in the desired range (0.17-0.55).")
    else:
        reasons.append(f"O/C ratio ({o_c_ratio:.2f}) is not in the desired range (0.17-0.55).")
    
    # Criterion 5: Furan ring check via SMARTS pattern.
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    has_furan = mol.HasSubstructMatch(furan_pattern)
    if has_furan:
        total_score += 1.0
        reasons.append("Furan ring found in the structure.")
    else:
        reasons.append("No furan ring found in the structure.")
    
    # Additional penalty: Check for a sugar substructure (pyranose ring) that might indicate glycosylation.
    # A simple SMARTS pattern for a pyranose ring:
    sugar_smarts = "[C@H]1OC(CO)C(O)C(O)C1O"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if mol.HasSubstructMatch(sugar_pattern):
        total_score -= 1.5
        reasons.append("Sugar moiety detected; subtracting 1.5 points as glycosides are typically not limonoids.")
    
    # Decide thresholds that differ if a furan ring is present or not.
    # Genuine limonoids with a furan ring: threshold = 4.0;
    # de‐furano limonoids (no furan bonus): threshold = 3.5.
    threshold = 4.0 if has_furan else 3.5
    
    overall = ""
    if total_score >= threshold:
        overall = f"Meets heuristic criteria for a limonoid (score {total_score:.1f}/~{threshold} threshold)."
        decision = True
    else:
        overall = f"Does NOT meet heuristic criteria for a limonoid (score {total_score:.1f}/~{threshold})."
        decision = False
    
    # Combine reasons
    reason_detail = overall + " Details: " + " ".join(reasons)
    return decision, reason_detail

# Example testing (when running as a script):
if __name__ == "__main__":
    # Test using Deacetylnomilin (a known limonoid with a furan ring)
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_limonoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)