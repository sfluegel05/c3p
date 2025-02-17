"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid class – tetranortriterpenoids that are highly oxygenated and either contain or are derived from
a 4,4,8-trimethyl-17-furanylsteroid skeleton. 
Heuristic criteria (each can contribute to a total score):
  1. Molecular weight between 350 and 650 Da.
  2. Polycyclic structure: at least 4 rings.
  3. Carbon count between 23 and 40.
  4. Oxygenation: oxygen/carbon ratio. Give 1 point if ≥ 0.18, or 0.5 if between 0.17 and 0.18.
  5. Furan ring presence (SMARTS "c1occc1").

A total score of 3.5 or greater classifies the molecule as a limonoid.
Note: This heuristic filter is not perfect—the thresholds were chosen to try to improve performance
over the previous attempt.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    
    Heuristic criteria used (each criterion contributing to a score):
      1. Molecular weight between 350 and 650 Da.
      2. Contains at least 4 rings.
      3. Carbon count between 23 and 40.
      4. High oxygenation: if O/C ratio >= 0.18 add 1 point; if in [0.17, 0.18) add 0.5.
      5. Contains a furan ring (SMARTS "c1occc1").
      
    If the total score is 3.5 or more (out of a maximum of 5), the molecule is classified as a limonoid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a limonoid, False otherwise.
        str: Detailed reasoning for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    reasons = []
    score = 0.0  # accumulate evidence for limonoid
    
    # Criterion 1: Molecular Weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if 350 <= mol_wt <= 650:
        score += 1
        reasons.append(f"MW {mol_wt:.1f} Da is within expected range (350-650).")
    else:
        reasons.append(f"MW {mol_wt:.1f} Da is outside the expected range (350-650).")
    
    # Criterion 2: Ring count (polycyclicity)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings >= 4:
        score += 1
        reasons.append(f"Number of rings ({num_rings}) is sufficient (>=4).")
    else:
        reasons.append(f"Only {num_rings} rings found (expected at least 4).")
    
    # Criterion 3: Carbon count check
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if 23 <= carbon_count <= 40:
        score += 1
        reasons.append(f"Carbon count ({carbon_count}) is within the expected range (23-40).")
    else:
        reasons.append(f"Carbon count ({carbon_count}) is outside the expected range (23-40).")
    
    # Criterion 4: Oxygenation (O/C ratio)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    oxygen_ratio = oxygen_count / carbon_count if carbon_count > 0 else 0
    if oxygen_ratio >= 0.18:
        score += 1
        reasons.append(f"O/C ratio ({oxygen_ratio:.2f}) is high (>=0.18).")
    elif 0.17 <= oxygen_ratio < 0.18:
        score += 0.5
        reasons.append(f"O/C ratio ({oxygen_ratio:.2f}) is borderline (>=0.17 but < 0.18); awarding half point.")
    else:
        reasons.append(f"O/C ratio ({oxygen_ratio:.2f}) is too low (expected >=0.17).")
    
    # Criterion 5: Furan ring check using SMARTS pattern for an aromatic furan
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    if mol.HasSubstructMatch(furan_pattern):
        score += 1
        reasons.append("Furan ring found in the structure.")
    else:
        reasons.append("No furan ring found.")
    
    # Final decision
    if score >= 3.5:
        overall = f"Meets heuristic criteria for a limonoid (score {score:.1f}/5)."
        return True, overall + " Details: " + " ".join(reasons)
    else:
        overall = f"Does NOT meet heuristic criteria for a limonoid (score {score:.1f}/5)."
        return False, overall + " Details: " + " ".join(reasons)

# Example testing when running as a script:
if __name__ == "__main__":
    # Example: using Deacetylnomilin (a limonoid example)
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_limonoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)