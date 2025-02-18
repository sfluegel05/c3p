"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
Definition: A 1-O-acylglycerophosphoethanolamine having (R)-configuration.
This program checks for:
  (1) A phosphoethanolamine headgroup (presence of an OCCN fragment where N is neutral),
  (2) Exactly one acyl ester (an ester “C(=O)O[C]” where the oxygen is not part of a phosphate ester),
  (3) A chiral center (with (R) configuration) that is part of the glycerol backbone (i.e. attached to a phosphorus).
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine (1-acyl-GPE)
    based on its SMILES string.

    The molecule must have:
      - A phosphoethanolamine headgroup (an OCCN fragment with a neutral nitrogen).
      - An acyl ester linkage connecting a fatty acid at sn-1 (exactly one ester of the form C(=O)O[C] that is not a phosphate ester).
      - At least one glycerol-associated chiral center with (R) configuration (tested as a chiral center bonded at least to one phosphorus).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 1‐acyl‐sn‐glycero‐3‐phosphoethanolamine, False otherwise.
        str: Explanation of the reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Molecule could not be sanitized: " + str(e)
    # assign stereochemistry if needed
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # 1. Check for phosphoethanolamine headgroup.
    # We look for an OCCN fragment (oxygen-carbon-carbon-nitrogen).
    headgroup_smarts = "OCCN"  
    headgroup_pat = Chem.MolFromSmarts(headgroup_smarts)
    headgroup_matches = mol.GetSubstructMatches(headgroup_pat)
    if not headgroup_matches:
        return False, "Phosphoethanolamine headgroup (OCCN) not found"
    # Verify that at least one matching nitrogen is neutral (atomic num 7, formal charge 0)
    headgroup_found = False
    for match in headgroup_matches:
        # match is a tuple of atom indices corresponding to O, C, C, N in order.
        n_atom = mol.GetAtomWithIdx(match[3])
        if n_atom.GetAtomicNum() == 7 and n_atom.GetFormalCharge() == 0:
            headgroup_found = True
            break
    if not headgroup_found:
        return False, "Phosphoethanolamine headgroup found, but nitrogen is not neutral"
    
    # 2. Check for the acyl ester linkage.
    # We use a SMARTS for an ester: a carbonyl (C(=O)) attached to an oxygen that's bonded to a carbon.
    # (Do not confuse phosphate esters, so we filter out matches where the O-attached carbon is bonded to P.)
    ester_smarts = "C(=O)O[C]"
    ester_pat = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pat)
    acyl_ester_count = 0
    for match in ester_matches:
        # match indices: 0->carbonyl C, 1->carbonyl O, 2->carbon attached to O.
        o_attached_idx = match[2]
        o_attached_atom = mol.GetAtomWithIdx(o_attached_idx)
        # If any neighbor of this atom is phosphorus, then this ester is likely part of a phosphate ester, so skip.
        is_phosphate_ester = any(neigh.GetAtomicNum() == 15 for neigh in o_attached_atom.GetNeighbors())
        if not is_phosphate_ester:
            acyl_ester_count += 1
    if acyl_ester_count != 1:
        return False, f"Found {acyl_ester_count} acyl ester linkage(s) (expected exactly 1 for 1-acyl lipid)"
    
    # 3. Check for a (R)-configured glycerol chiral center.
    # We expect the glycerol backbone to contain a chiral center (typically the sn-2 carbon)
    # that is bonded to the phosphate. So we search for chiral centers, then check if any (R) center
    # has at least one neighboring phosphorus.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    r_center_found = False
    for idx, config in chiral_centers:
        # Look at neighbors to see if one is a phosphorus atom.
        atom = mol.GetAtomWithIdx(idx)
        neighbors = atom.GetNeighbors()
        bonded_to_P = any(neigh.GetAtomicNum() == 15 for neigh in neighbors)
        if bonded_to_P and config.upper() == "R":
            r_center_found = True
            break
    if not r_center_found:
        return False, "No (R)-configured chiral center (bonded to phosphorus) found; required glycerol configuration missing"
    
    return True, "Molecule passes all tests: it has a phosphoethanolamine headgroup, exactly one acyl ester linkage (sn-1), and an (R)-configured glycerol chiral center"

# Optional testing block – can be removed if not desired.
if __name__ == "__main__":
    # Test a valid molecule: 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
    test_smiles_valid = "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"
    result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(test_smiles_valid)
    print("Test valid molecule:", result)
    print("Reason:", reason)
    
    # Example of a false positive: a diacyl phospholipid (e.g., phosphatidylcholine) that should not match.
    test_smiles_false = "C([C@@](COC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)=O)[H])OP([O-])(=O)OCC[N+](C)(C)C"
    result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(test_smiles_false)
    print("Test false positive molecule:", result)
    print("Reason:", reason)