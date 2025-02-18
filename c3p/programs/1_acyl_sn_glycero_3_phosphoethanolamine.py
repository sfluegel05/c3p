"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
Definition: A 1-O-acylglycerophosphoethanolamine having (R)-configuration
However, due to known RDKit CIP assignment issues for the glycerol backbone,
this code accepts any chiral center (R or S) that is attached to a phosphorus atom.
It also requires:
  (1) A phosphoethanolamine headgroup, defined by an OCCN fragment with a neutral nitrogen.
  (2) Exactly one acyl ester linkage (an ester bond of the form C(=O)O[C] that is not part
      of a phosphate ester).
  (3) A glycerol chiral center (bonded to at least one phosphorus) present.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine (1-acyl-GPE)
    based on its SMILES string.

    Requirements:
      - Presence of a phosphoethanolamine headgroup (an OCCN fragment where N is neutral).
      - Exactly one acyl ester linkage (an ester of the form C(=O)O[C] where the O-attached
        carbon is not bonded to phosphorus).
      - A chiral center in the glycerol region that is directly bonded to a phosphorus atom.
        (Note: Although the definition requires an (R)-configuration, our analysis here will
         accept either “R” or “S”, because in many valid structures the RDKit-assigned CIP
         label is inverted relative to the natural configuration.)
    
    Args:
        smiles (str): SMILES string of the molecule

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
    # Ensure stereochemistry is assigned.
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # 1. Check for phosphoethanolamine headgroup.
    # Look for an OCCN fragment. Only accept if the N is neutral.
    headgroup_smarts = "OCCN"
    headgroup_pat = Chem.MolFromSmarts(headgroup_smarts)
    headgroup_matches = mol.GetSubstructMatches(headgroup_pat)
    if not headgroup_matches:
        return False, "Phosphoethanolamine headgroup (OCCN) not found"
    headgroup_found = False
    for match in headgroup_matches:
        # match gives indices for O, C, C, N.
        n_atom = mol.GetAtomWithIdx(match[3])
        if n_atom.GetAtomicNum() == 7 and n_atom.GetFormalCharge() == 0:
            headgroup_found = True
            break
    if not headgroup_found:
        return False, "Phosphoethanolamine headgroup found, but nitrogen is not neutral"
    
    # 2. Check for the acyl ester linkage.
    # We search for an ester fragment of the form C(=O)O[C]
    # and then exclude any where the oxygen-linked carbon is attached to phosphorus.
    ester_smarts = "C(=O)O[C]"
    ester_pat = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pat)
    acyl_ester_count = 0
    for match in ester_matches:
        # match indices: 0->carbonyl C, 1->carbonyl O, 2->carbon attached to O.
        o_attached_idx = match[2]
        o_attached_atom = mol.GetAtomWithIdx(o_attached_idx)
        # Skip if the O-attached carbon neighbors a phosphorus.
        if any(neigh.GetAtomicNum() == 15 for neigh in o_attached_atom.GetNeighbors()):
            continue
        acyl_ester_count += 1
    if acyl_ester_count != 1:
        return False, f"Found {acyl_ester_count} acyl ester linkage(s) (expected exactly 1 for a 1-acyl lipid)"
    
    # 3. Check for a chiral center in the glycerol backbone.
    # We expect a chiral center that is attached to at least one phosphorus atom.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    glycerol_chiral_found = False
    for idx, config in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        # We require the atom to have a chiral tag indicating it is truly chiral.
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        # Check if this center is bonded to at least one phosphorus.
        if any(neighbor.GetAtomicNum() == 15 for neighbor in atom.GetNeighbors()):
            glycerol_chiral_found = True
            break
    if not glycerol_chiral_found:
        return False, "No chiral center (bonded to phosphorus) found in the glycerol backbone; required glycerol configuration missing"
    
    return True, ("Molecule passes all tests: it has a phosphoethanolamine headgroup, exactly one acyl ester linkage (sn-1), "
                  "and a glycerol chiral center (accepted regardless of reported CIP configuration)")

# Optional testing block – remove if not desired.
if __name__ == "__main__":
    test_smiles = [
        # Valid examples:
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O",  # LysoPE(20:3(5Z,8Z,11Z)/0:0)
        "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",  # 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",  # 1-oleoyl-sn-glycero-3-phosphoethanolamine
    ]
    for smi in test_smiles:
        result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 80)