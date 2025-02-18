"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid
Definition: Any member of the class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage
to a carbohydrate part (usually a mono-, di- or tri-saccharide). Some bacterial glycolipids have the sugar part acylated,
and the glycerol component may be absent.
This program uses several heuristic tests:
  1. It looks for a sugar moiety by scanning for a nonaromatic ring of size 5 or 6 that contains exactly one oxygen atom.
  2. It looks for a glycosidic (ether) linkage connecting a non‐ring carbon to a ring carbon.
  3. It looks for an acyl linkage (ester or amide) that could attach a fatty acid.
  4. It checks for the presence of a long aliphatic chain (heuristically, 8 or more consecutive carbon atoms).
  5. It checks that the molecular weight is above 500 Da and that there are a minimum number of rotatable bonds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def has_sugar(mol):
    """
    Looks for a candidate sugar ring.
    We define a sugar ring as a nonaromatic ring of size 5 or 6 having exactly one oxygen atom in the ring.
    (Many carbohydrate rings are pyranoses/furanoses with one ring oxygen.)
    """
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) not in [5, 6]:
            continue
        # Skip rings that are aromatic (sugars are not aromatic)
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Count oxygen atoms in the ring 
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxy_count == 1:
            return True
    return False

def has_long_aliphatic_chain(smiles):
    """
    Heuristic: check if the SMILES string contains 8 or more consecutive 'C' characters.
    To avoid being misled by bond characters, we remove backslash and forward slash.
    """
    clean_smiles = smiles.replace("\\", "").replace("/", "")
    return "CCCCCCCC" in clean_smiles

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    The algorithm uses the following heuristic tests:
      1. Identify a sugar moiety via a nonaromatic five-/six-membered ring with exactly one oxygen.
      2. Confirm a glycosidic linkage exists by finding an ether-bond [C!R]-O-[C R] (a non-ring to ring connection).
      3. Identify at least one acyl linkage (ester or amide) that can anchor a fatty acyl chain.
      4. Verify that there is a long aliphatic chain (at least 8 consecutive carbon characters).
      5. Ensure that molecular weight (>500 Da) and the number of rotatable bonds (>= 3) are in an acceptable range.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the first element indicates if the molecule is classified as a glycolipid,
                     and the second element gives the reason for the classification or failure.
    """
    # Parse the SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Test 1: Check for a sugar ring moiety
    if not has_sugar(mol):
        return False, "No carbohydrate (sugar ring) moiety found"
    
    # Test 2: Check for a glycosidic (ether) linkage connecting a non-ring carbon to a ring carbon.
    # This SMARTS pattern looks for a bond O linking a non-ring carbon ([C;!R]) and a ring carbon ([C;R]).
    glyco_pattern = Chem.MolFromSmarts("[C;!R]-O-[C;R]")
    if glyco_pattern is None or not mol.HasSubstructMatch(glyco_pattern):
        return False, "No glycosidic linkage (ether bond between a non‐ring and a ring atom) found"
    
    # Test 3: Check for acyl linkages.
    # Ester linkage: [OX2][CX3](=O)[#6]
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)[#6]")
    if ester_pattern is None:
        return False, "Error in ester SMARTS pattern"
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Amide linkage: NC(=O)[#6]
    amide_pattern = Chem.MolFromSmarts("NC(=O)[#6]")
    if amide_pattern is None:
        return False, "Error in amide SMARTS pattern"
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    if (len(ester_matches) + len(amide_matches)) < 1:
        return False, "No acyl linkage (ester or amide) found that could anchor a fatty acid chain"
    
    # Test 4: Check for the presence of a long fatty acyl chain
    if not has_long_aliphatic_chain(smiles):
        return False, "No long aliphatic (fatty acid) chain detected (need at least 8 consecutive carbons)"
    
    # Test 5: Check molecular weight (typically >500 Da for glycolipids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a glycolipid"
    
    # Test 6: Check that there are enough rotatable bonds (>=3)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Not enough rotatable bonds to support the presence of long fatty acyl chains"
    
    return True, "Molecule contains a sugar ring linked via a glycosidic bond to an acylated lipid chain consistent with a glycolipid structure"

# For testing purposes, you may uncomment the following:
# if __name__ == "__main__":
#     test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2OCCCCCCCC\\C=C/CCCCCCCC(O)=O"
#     result, reason = is_glycolipid(test_smiles)
#     print("Is glycolipid?", result)
#     print("Reason:", reason)