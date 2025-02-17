"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid class – tetranortriterpenoids that are highly oxygenated and either contain 
or are derived from a precursor with a 4,4,8-trimethyl-17-furanylsteroid skeleton.
Heuristic criteria (each can contribute to a total score):
  1. Molecular weight between 350 and 750 Da. (1 point)
  2. Polycyclic structure:
       - At least 4 rings yields 1 point.
       - If there are 6 or more rings, add a further bonus of 0.5.
  3. Carbon count between 23 and 40. (1 point)
  4. Oxygenation: If O/C ratio is ≥ 0.17, add 1 point.
  5. Furan ring presence (SMARTS "c1occc1"). (1 point)
  
A total score of 4.0 or more (out of ~5.0–5.5 maximum) is taken as evidence that the molecule
belongs to the limonoid class.
Note:
  – Some limonoids (for example, defurano compounds) might lose the furan bonus.
  – This heuristic is an imperfect filter, but our modifications (wider MW range, ring bonus, and
    a full point for oxygenation) were intended to improve performance.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines whether a molecule can be considered a limonoid based on heuristic criteria.
    
    The criteria (with associated points) are:
      1. Molecular Weight between 350 and 750 Da (1 point).
      2. Number of rings: ≥4 (1 point) and if ≥6 rings add an extra bonus of 0.5.
      3. Carbon count between 23 and 40 (1 point).
      4. Oxygenation: O/C ratio ≥ 0.17 adds 1 point.
      5. Presence of a furan ring (SMARTS "c1occc1") adds 1 point.
      
    A total score ≥ 4.0 is taken as evidence toward limonoid classification.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a limonoid, False otherwise.
        str: Reason detailing how the score was computed and which criteria were met.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    reasons = []
    total_score = 0.0  # accumulate heuristic score

    # Criterion 1: Molecular weight between 350 and 750 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if 350 <= mol_wt <= 750:
        total_score += 1.0
        reasons.append(f"MW {mol_wt:.1f} Da is within expected range (350-750).")
    else:
        reasons.append(f"MW {mol_wt:.1f} Da is outside the expected range (350-750).")
    
    # Criterion 2: Count rings (polycyclic).
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings >= 4:
        total_score += 1.0
        reasons.append(f"Number of rings ({num_rings}) is sufficient (>=4).")
        if num_rings >= 6:
            total_score += 0.5
            reasons.append("Bonus: Rings count is high (>=6), adding extra 0.5.")
    else:
        reasons.append(f"Only {num_rings} rings found (expected at least 4).")
    
    # Criterion 3: Carbon count.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if 23 <= carbon_count <= 40:
        total_score += 1.0
        reasons.append(f"Carbon count ({carbon_count}) is within expected range (23-40).")
    else:
        reasons.append(f"Carbon count ({carbon_count}) is outside expected range (23-40).")
    
    # Criterion 4: Oxygenation level (Oxygen-to-Carbon ratio).
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    o_c_ratio = oxygen_count / carbon_count if carbon_count > 0 else 0
    if o_c_ratio >= 0.17:
        total_score += 1.0
        reasons.append(f"O/C ratio ({o_c_ratio:.2f}) is high (>=0.17).")
    else:
        reasons.append(f"O/C ratio ({o_c_ratio:.2f}) is too low (expected ≥0.17).")
    
    # Criterion 5: Furan ring check.
    # SMARTS for an aromatic furan ring.
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    if mol.HasSubstructMatch(furan_pattern):
        total_score += 1.0
        reasons.append("Furan ring found in the structure.")
    else:
        reasons.append("No furan ring found in the structure.")
    
    # Final decision.
    # (A total score of 4.0 or more is taken as evidence toward a limonoid.)
    overall = ""
    if total_score >= 4.0:
        overall = f"Meets heuristic criteria for a limonoid (score {total_score:.1f}/~5.0)."
        return True, overall + " Details: " + " ".join(reasons)
    else:
        overall = f"Does NOT meet heuristic criteria for a limonoid (score {total_score:.1f}/~5.0)."
        return False, overall + " Details: " + " ".join(reasons)


# Example testing when running as a script:
if __name__ == "__main__":
    # Test using Deacetylnomilin (a known limonoid)
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_limonoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)