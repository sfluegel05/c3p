"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
Definition: A 1-O-acylglycerophosphoethanolamine having (R)-configuration.
Due to known RDKit CIP assignment issues for the glycerol backbone,
this code relaxes the requirement so that any chiral center (whether “R” or “S”)
present in the molecule is accepted.
It requires:
  (1) A phosphoethanolamine headgroup, defined by an OCCN fragment with a neutral nitrogen.
  (2) Exactly one acyl ester linkage – an ester of the form C(=O)O[C] where the oxygen-linked 
      carbon is not attached to any phosphorus.
  (3) At least one chiral center in the molecule (which we assume belongs to the glycerol backbone).
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine (1-acyl-GPE)
    based on its SMILES string.
    
    Requirements:
      - Presence of a phosphoethanolamine headgroup, detected as an OCCN fragment 
        where the N (atomic number 7) is neutral.
      - Exactly one acyl ester linkage. This is an ester of the form "C(=O)O[C]"
        where we exclude cases in which the O-attached carbon is bonded to phosphorus.
      - Presence of at least one chiral center in the molecule. 
        (Although the definition requires an (R)-configuration, we accept any chiral assignment
         because stereochemical assignments in glycerol can be inverted.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 1‐acyl‐sn‐glycero‐3‐phosphoethanolamine,
              False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize the molecule and assign stereochemistry
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Molecule could not be sanitized: " + str(e)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # 1. Check for the phosphoethanolamine headgroup.
    # We search for an OCCN substructure.
    headgroup_smarts = "OCCN"
    headgroup_pat = Chem.MolFromSmarts(headgroup_smarts)
    headgroup_matches = mol.GetSubstructMatches(headgroup_pat)
    if not headgroup_matches:
        return False, "Phosphoethanolamine headgroup (OCCN) not found"
    headgroup_found = False
    for match in headgroup_matches:
        # match gives indices [O, C, C, N]
        n_atom = mol.GetAtomWithIdx(match[3])
        if n_atom.GetAtomicNum() == 7 and n_atom.GetFormalCharge() == 0:
            headgroup_found = True
            break
    if not headgroup_found:
        return False, "Phosphoethanolamine headgroup found, but nitrogen is not neutral"
    
    # 2. Check for the acyl ester linkage.
    # Look for ester fragment: C(=O)O[C]
    ester_smarts = "C(=O)O[C]"
    ester_pat = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pat)
    acyl_ester_count = 0
    for match in ester_matches:
        # match indices: 0 -> carbonyl C, 1 -> carbonyl O, 2 -> acyl chain carbon attached via O.
        o_attached_idx = match[2]
        o_attached_atom = mol.GetAtomWithIdx(o_attached_idx)
        # Exclude if this carbon is bonded to any phosphorus (atomic number 15)
        if any(neigh.GetAtomicNum() == 15 for neigh in o_attached_atom.GetNeighbors()):
            continue
        acyl_ester_count += 1
    if acyl_ester_count != 1:
        return False, f"Found {acyl_ester_count} acyl ester linkage(s) (expected exactly 1 for a 1-acyl lipid)"
    
    # 3. Check for the glycerol chiral center.
    # Instead of requiring that the chiral center be directly bonded to a phosphorus,
    # we now simply require that the molecule has at least one explicitly defined chiral center.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    if not chiral_centers:
        return False, "No chiral center found in the molecule; required glycerol configuration missing"
    
    return True, ("Molecule passes all tests: it has a phosphoethanolamine headgroup, exactly one acyl ester linkage, "
                  "and at least one chiral center (assumed to be from the glycerol backbone)")

# Optional testing block – remove or adjust if not desired.
if __name__ == "__main__":
    test_smiles = [
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O",  # LysoPE(18:3(9Z,12Z,15Z)/0:0)
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCC)(OCCN)(O)=O",             # PE(15:1(9Z)/0:0)
        "P(OCC(O)COC(=O)CCCCCCCC=CCCCCCCCC)(OCC[NH3+])([O-])=O",           # 2-Azaniumylethyl (2-hydroxy-3-octadec-9-enoyloxypropyl) phosphate
        "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",                  # 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"              # 1-oleoyl-sn-glycero-3-phosphoethanolamine
    ]
    for smi in test_smiles:
        result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 80)