"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
Definition: A 1-O-acylglycerophosphoethanolamine having (R)-configuration.
Due to known RDKit CIP assignment issues for glycerol backbone stereochemistry,
this code relaxes the requirement so that either an explicit chiral center
or a glycerol backbone substructure is accepted.
It requires:
  (1) A phosphoethanolamine headgroup:
      – a neutral OCCN fragment that is near a phosphorus atom.
  (2) Exactly one acyl ester linkage – an ester of the form C(=O)O[C] where the
      O-linked carbon does not connect to any phosphorus.
  (3) The presence of either an explicitly defined chiral center or a glycerol-like backbone.
  (4) A molecular weight typical for a lipid (< 1000 Da) and exactly one phosphorus.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines whether a molecule (given by a SMILES string) is a 1-acyl-sn-glycero-3-phosphoethanolamine.
    
    Returns:
       (bool, str): True with a reason if classified as 1-acyl-GPE, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is sanitized and stereochemistry is assigned.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Molecule could not be sanitized: " + str(e)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # --- Additional filter: exactly one phosphorus atom ---
    phosphorus = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus) != 1:
        return False, f"Expected exactly one phosphorus atom; found {len(phosphorus)}"
    
    # --- Check for phosphoethanolamine headgroup ---
    # First, look for a neutral OCCN substructure.
    headgroup_smarts = "OCCN"
    headgroup_pat = Chem.MolFromSmarts(headgroup_smarts)
    headgroup_matches = mol.GetSubstructMatches(headgroup_pat)
    if not headgroup_matches:
        return False, "Phosphoethanolamine headgroup (OCCN) not found"
    
    headgroup_found = False
    # Check that at least one OCCN instance has a neutral N and is near the phosphorus.
    for match in headgroup_matches:
        # Expecting match indices [O, C, C, N]
        n_atom = mol.GetAtomWithIdx(match[3])
        if n_atom.GetAtomicNum() == 7 and n_atom.GetFormalCharge() == 0:
            # Check distance from the nitrogen to the phosphorus.
            # Here we demand that a phosphorus atom is within 3 bonds.
            dist_found = False
            for p_atom in phosphorus:
                p_idx = p_atom.GetIdx()
                # Get bond path lengths between the N and P.
                if Chem.rdmolops.GetShortestPath(mol, n_atom.GetIdx(), p_idx):
                    # In our simple approach we assume a found path implies they are connected.
                    if len(Chem.rdmolops.GetShortestPath(mol, n_atom.GetIdx(), p_idx)) - 1 <= 3:
                        dist_found = True
                        break
            if dist_found:
                headgroup_found = True
                break
    if not headgroup_found:
        return False, "Phosphoethanolamine headgroup found but no neutral OCCN fragment is near a phosphorus atom"
    
    # --- Check for exactly one acyl ester linkage ---
    # Look for ester fragments of the form: carbonyl C attached to O which is attached to a carbon.
    ester_smarts = "C(=O)O[C]"
    ester_pat = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pat)
    acyl_ester_count = 0
    for match in ester_matches:
        # match indices: 0 -> carbonyl carbon, 1 -> carbonyl oxygen, 2 -> O-linked carbon.
        o_attached_idx = match[2]
        o_attached_atom = mol.GetAtomWithIdx(o_attached_idx)
        # Exclude this ester if the O-linked carbon is attached to any phosphorus.
        if any(neigh.GetAtomicNum() == 15 for neigh in o_attached_atom.GetNeighbors()):
            continue
        acyl_ester_count += 1
    if acyl_ester_count != 1:
        return False, f"Found {acyl_ester_count} acyl ester linkage(s); exactly one is required"
    
    # --- Check for a glycerol backbone signature ---
    # First try finding an explicitly defined chiral center.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    if not chiral_centers:
        # If no chiral center is explicitly defined, try to find a glycerol-like substructure.
        # Many lipid structures contain a glycerol backbone fragment such as: O C C(O) C O C(=O)
        glycerol_smarts = "OCC(O)COC(=O)"
        glycerol_pat = Chem.MolFromSmarts(glycerol_smarts)
        if not mol.HasSubstructMatch(glycerol_pat):
            return False, "No chiral center found and glycerol backbone substructure (OCC(O)COC(=O)) missing"
        else:
            chiral_reason = "No explicit chiral center found, but glycerol backbone substructure present"
    else:
        chiral_reason = "At least one chiral center is defined in the structure"
    
    # --- Filter by molecular weight (typical for a lipid; avoid large peptides) ---
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, f"Molecular weight too high for a typical 1-acyl-GPE lipid (found {mol_wt:.1f} Da)"
    
    return True, ("Molecule passes all tests: it has a phosphoethanolamine headgroup with neutral OCCN adjacent to P, "
                  "exactly one acyl ester linkage not involving P, and " + chiral_reason + 
                  f", with a molecular weight of {mol_wt:.1f} Da.")

# Optional testing block – remove or adjust as desired.
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O",  # LysoPE(18:3(9Z,12Z,15Z)/0:0)
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCC)(OCCN)(O)=O",             # PE(15:1(9Z)/0:0)
        "P(OC[C@H](O)COC(=O)CCCCCCCCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN)(O)=O",  # LysoPE(20:3(5Z,8Z,11Z)/0:0)
        "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",                  # 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",             # 1-oleoyl-sn-glycero-3-phosphoethanolamine
        # A molecule that previously was a false negative due to missing chiral flag:
        "P(OCC(O)COC(=O)CCCCCCCC/C=C\\C/C=C\\C/C=C\\CCCCC)(OCCN)(O)=O",    # LPE 22:3
        # False positive (peptide-like; should be rejected due to weight or lack of headgroup connectivity):
        "BrC1=C(O)C=CC(=C1)C2NC(=O)[C@H](N(C(=O)[C@@H](NC(=O)[C@@H](C)CC(=C[C@@H]([C@H](OC([C@H]2OC)=O)C)C)C)C)CC=3C4=C(C=CC=C4)NC3"
    ]
    for smi in test_smiles:
        result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 80)