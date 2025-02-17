"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid
Definition: Any member of the class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage to a carbohydrate part (usually a mono-, di- or tri-saccharide). Some substances (e.g. bacterial glycolipids) may have the sugar part acylated by one or more fatty acids and the glycerol part may be absent.
This program uses several heuristic tests:
  1. It looks for a sugar ring by looping over all rings (of size 5–6) and checking for predominance of carbon and oxygen.
  2. It looks for at least one acyl connection (ester or amide).
  3. It checks if a long lipid chain is present (heuristically, a substring of eight consecutive carbon atoms in the SMILES).
  4. It checks that the molecular weight is above 500 Da and there are enough rotatable bonds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def has_sugar(mol):
    """
    Checks for presence of a sugar-like ring. We define a candidate sugar ring as a ring of size 5 or 6 where
    at least half the atoms are oxygen or carbon and at least two atoms are oxygen.
    """
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) not in [5, 6]:
            continue
        oxy_count = 0
        valid = True
        for idx in ring:
            atomic_num = mol.GetAtomWithIdx(idx).GetAtomicNum()
            if atomic_num not in (6, 8):  # only carbon and oxygen allowed
                valid = False
                break
            if atomic_num == 8:
                oxy_count += 1
        # require at least two oxygens and ring is not only carbons
        if valid and oxy_count >= 2:
            return True
    return False

def has_long_aliphatic_chain(smiles):
    """
    A simple heuristic: check if the SMILES string contains 8 or more consecutive 'C' characters.
    This is not mechanistically perfect but catches many long fatty acid chains.
    """
    # remove bond characters to help the string search:
    clean_smiles = smiles.replace("\\", "").replace("/", "")
    if "CCCCCCCC" in clean_smiles:
        return True
    return False

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    The algorithm uses several heuristic tests:
      1. It checks for a sugar-like moiety by scanning for rings (5–6 members) rich in O atoms.
      2. It checks for an acyl bond (ester or amide link) that could attach a fatty acid.
      3. It verifies the molecule has a long aliphatic chain (by checking for 8+ consecutive carbons).
      4. It checks that the molecular weight is moderately high (e.g. >500 Da).
      5. It also makes sure there are enough rotatable bonds (>= 3) to allow long chain flexibility.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a glycolipid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for sugar moiety using our helper function.
    if not has_sugar(mol):
        return False, "No carbohydrate (sugar ring) moiety found"
    
    # Check for acyl linkage via ester pattern.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)[#6]")
    if ester_pattern is None:
        return False, "Error in ester SMARTS pattern"
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check for amide pattern.
    amide_pattern = Chem.MolFromSmarts("NC(=O)[#6]")
    if amide_pattern is None:
        return False, "Error in amide SMARTS pattern"
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    if (len(ester_matches) + len(amide_matches)) < 1:
        return False, "No acyl linkage (ester or amide) found that could anchor a fatty acid chain"
    
    # Check for the presence of a long fatty acyl chain using a simple string heuristic.
    if not has_long_aliphatic_chain(smiles):
        return False, "No long aliphatic (fatty acid) chain detected (need at least 8 consecutive carbons)"
    
    # Check molecular weight - glycolipids are often >500 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a glycolipid"
    
    # Check if we have a reasonable number of rotatable bonds.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Not enough rotatable bonds to support presence of long fatty acyl chains"
    
    return True, "Molecule contains a sugar moiety, acyl linkage(s), and a long fatty acid chain consistent with a glycolipid structure"

# For testing purposes (uncomment to run tests):
# if __name__ == "__main__":
#     test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2OCCCCCCCC\\C=C/CCCCCCCC(O)=O"
#     result, reason = is_glycolipid(test_smiles)
#     print("Is glycolipid?", result)
#     print("Reason:", reason)