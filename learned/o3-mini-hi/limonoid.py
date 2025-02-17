"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid class – tetranortriterpenoids that are highly oxygenated and either contain or are derived from
a 4,4,8-trimethyl-17-furanylsteroid skeleton. This is a heuristic filter using multiple criteria.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    
    Heuristics used (each satisfied criterion adds 1 point):
      1. Molecular weight in the range 350-600 Da.
      2. Polycyclic structure (at least 4 rings).
      3. Carbon count between 23 and 40.
      4. High oxygenation (oxygen/carbon ratio >= 0.20).
      5. Contains a furan ring (SMARTS: "c1occc1").
         Note: Defurano derivatives may be accepted if the other criteria are strong.
      
    If the total score is 4 or more out of 5, the molecule is classified as a limonoid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the heuristic criteria for a limonoid, False otherwise.
        str: Detailed reasoning for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    reasons = []
    score = 0  # accumulate evidence for limonoid
    
    # 1. Molecular Weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if 350 <= mol_wt <= 600:
        score += 1
        reasons.append(f"MW {mol_wt:.1f} Da is within expected range (350-600).")
    else:
        reasons.append(f"MW {mol_wt:.1f} Da is outside the expected range (350-600).")
    
    # 2. Polycyclic (ring count) check
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings >= 4:
        score += 1
        reasons.append(f"Number of rings ({num_rings}) is sufficient (>=4).")
    else:
        reasons.append(f"Only {num_rings} rings found (expected at least 4).")
    
    # 3. Carbon count check
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if 23 <= carbon_count <= 40:
        score += 1
        reasons.append(f"Carbon count ({carbon_count}) is within the expected range (23-40).")
    else:
        reasons.append(f"Carbon count ({carbon_count}) is outside the expected range (23-40).")
    
    # 4. Oxygenation check: compute oxygen/carbon ratio
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    oxygen_ratio = oxygen_count / carbon_count if carbon_count > 0 else 0
    if oxygen_ratio >= 0.20:
        score += 1
        reasons.append(f"O/C ratio ({oxygen_ratio:.2f}) is high (>=0.20).")
    else:
        reasons.append(f"O/C ratio ({oxygen_ratio:.2f}) is too low (expected >= 0.20).")
    
    # 5. Furan ring check
    # SMARTS for a furan ring: aromatic 5-membered ring with an oxygen.
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    has_furan = mol.HasSubstructMatch(furan_pattern)
    if has_furan:
        score += 1
        reasons.append("Furan ring found in the structure.")
    else:
        reasons.append("No furan ring found.")
    
    # Final decision based on total score.
    if score >= 4:
        overall = "Meets heuristic criteria for a limonoid (score {}/5).".format(score)
        return True, overall + " Details: " + " ".join(reasons)
    else:
        overall = "Does NOT meet heuristic criteria for a limonoid (score {}/5).".format(score)
        return False, overall + " Details: " + " ".join(reasons)


# Example testing when running as a script:
if __name__ == "__main__":
    # Example using Deacetylnomilin (a limonoid example)
    test_smiles = "O1[C@@]23[C@]4([C@@]([C@@]5([C@@](CC4=O)(C(OC(=O)CC5O)(C)C)[H])C)(CC[C@]2([C@@H](OC(=O)[C@@]13[H])C=6C=COC6)C)[H])C"
    result, reason = is_limonoid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)