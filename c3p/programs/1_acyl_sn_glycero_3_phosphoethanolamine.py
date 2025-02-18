"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: CHEBI:1-acyl-sn-glycero-3-phosphoethanolamine
Definition: A 1-O-acylglycerophosphoethanolamine having (R)-configuration.
This classifier requires:
  (1) Exactly one phosphorus atom.
  (2) A phosphoethanolamine headgroup in which an "OCCN" fragment is present;
      the oxygen of that fragment must be directly bonded to the phosphorus.
  (3) Exactly one acyl ester linkage (i.e. a C(=O)O–C fragment) where the O-linked carbon
      is not linked to any phosphorus.
  (4) A glycerol backbone substructure as seen in many glycerophospholipids.
  (5) A molecular weight below 1000 Da.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdmolops

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines whether a molecule (given via its SMILES string) is a 1-acyl-sn-glycero-3-phosphoethanolamine.
    
    Returns:
       bool: True if the molecule meets the criteria, False otherwise.
       str: The reason for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanitize and assign stereochemistry.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Molecule could not be sanitized: " + str(e)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # (1) Exactly one phosphorus atom.
    phosphorus = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus) != 1:
        return False, f"Expected exactly one phosphorus atom; found {len(phosphorus)}"
    
    # (2) Phosphoethanolamine headgroup check.
    # Look for an OCCN fragment.
    headgroup_pat = Chem.MolFromSmarts("OCCN")
    headgroup_matches = mol.GetSubstructMatches(headgroup_pat)
    headgroup_found = False
    for match in headgroup_matches:
        # The match is a tuple of atom indices: (O, C, C, N)
        o_idx = match[0]  # the oxygen that should be connected to phosphorus
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Check that one neighbor of this oxygen is a phosphorus.
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 15:
                # Also check that the nitrogen is neutral.
                n_atom = mol.GetAtomWithIdx(match[3])
                if n_atom.GetAtomicNum() == 7 and n_atom.GetFormalCharge() == 0:
                    headgroup_found = True
                    break
        if headgroup_found:
            break
    if not headgroup_found:
        return False, "Phosphoethanolamine headgroup (neutral OCCN attached directly to P) not found"
    
    # (3) Exactly one acyl ester linkage.
    # Look for ester fragments of the form: carbonyl C attached to oxygen attached to another carbon.
    ester_pat = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pat)
    acyl_ester_count = 0
    for match in ester_matches:
        # Match indices: 0 -> carbonyl carbon; 1 -> carbonyl oxygen; 2 -> O-linked carbon.
        o_linked_atom = mol.GetAtomWithIdx(match[2])
        # Exclude if this O-linked carbon is attached to any phosphorus.
        if any(neigh.GetAtomicNum() == 15 for neigh in o_linked_atom.GetNeighbors()):
            continue
        acyl_ester_count += 1
    if acyl_ester_count != 1:
        return False, f"Found {acyl_ester_count} acyl ester linkage(s); exactly one is required"
    
    # (4) Detect glycerol backbone.
    # Many glycerophospholipids have a backbone fragment like: O-C-C(O)-C-O-P .
    # We use a SMARTS pattern "OCC(O)COP" (ignoring chirality markers) as a heuristic.
    glycerol_pat = Chem.MolFromSmarts("OCC(O)COP")
    if not mol.HasSubstructMatch(glycerol_pat):
        return False, "Glycerol backbone substructure (OCC(O)COP) not detected"
    
    # (5) Molecular weight typical for a lipid: below 1000 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, f"Molecular weight too high for a typical 1-acyl-GPE lipid (found {mol_wt:.1f} Da)"
    
    reason = (f"Molecule passes all tests: it contains a phosphoethanolamine headgroup (OCCN directly linked to P), "
              f"exactly one acyl ester linkage not bonded to P, a glycerol backbone substructure (OCC(O)COP), "
              f"and has a molecular weight of {mol_wt:.1f} Da.")
    return True, reason


# Optional testing block.
if __name__ == "__main__":
    # A list of example SMILES strings.
    test_smiles = [
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O",  # LysoPE(18:3(9Z,12Z,15Z)/0:0)
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCC)(OCCN)(O)=O",             # PE(15:1(9Z)/0:0)
        "P(OCC(O)COC(=O)CCCCCCCC=CCCCCCCCC)(OCC[NH3+])([O-])=O",           # 2-Azaniumylethyl (2-hydroxy-3-octadec-9-enoyloxypropyl) phosphate
        "P(OC[C@H](O)COC(=O)CCCCCCCCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN)(O)=O",  # LysoPE(20:3(5Z,8Z,11Z)/0:0)
        "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",                  # 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",             # 1-oleoyl-sn-glycero-3-phosphoethanolamine
        "C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(O)[H])OP(=O)(O)OCCN",         # PE(18:1(11Z)/0:0)
    ]
    for smi in test_smiles:
        result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 80)