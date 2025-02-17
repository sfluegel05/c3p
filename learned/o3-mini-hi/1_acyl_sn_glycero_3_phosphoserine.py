"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
#!/usr/bin/env python3
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine compounds
Defined as: An sn-glycerophosphoserine compound having an acyl substituent 
at the 1-hydroxy position.
The algorithm requires:
  (1) that the molecule contains a phosphoserine head group,
  (2) that a phosphorus atom is present,
  (3) that exactly one acyl ester branch from phosphorus is found that 
      represents the expected acyl substitution at the sn-1 position.
"""

from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines whether a molecule is a 1-acyl-sn-glycero-3-phosphoserine compound.
    The molecule must contain:
      (a) a phosphoserine head group (identified by a fragment like OCC(N)C(O)=O directly attached to phosphorus),
      (b) a phosphorus atom,
      (c) exactly one acyl ester substituent on the sn-1 (primary) hydroxyl branch.
    The acyl ester branch is assumed to follow the pattern:
      P - O - C - C(O) - C - O - C(=O) - R  (with R a carbon fragment)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple; True and explanation if the molecule belongs to the class,
                     otherwise False with a reason.
    """
    # Parse the SMILES string into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (1) Check that a phosphorus atom is present.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found"
    
    # (2) Check for the phosphoserine head group.
    # We require a branch directly attached to phosphorus matching a serine fragment.
    # This SMARTS uses the phosphorus atom at the “root” and expects one branch of form: 
    # O-C-C(N)C(O)=O (ignoring chirality specifications so variants with @ are allowed).
    serine_smarts = "P(OCC(N)C(O)=O)"
    serine_frag = Chem.MolFromSmarts(serine_smarts)
    if not mol.HasSubstructMatch(serine_frag):
        return False, "Phosphoserine head group not found"
    
    # (3) Check for the acyl ester substituent at the sn-1 hydroxyl.
    # We require that one substituent from phosphorus follows the expected pattern:
    #   P - O - C - C(O) - C - O - C(=O) - [#6]
    # Note: We ignore chirality so that variants with or without @ symbols match.
    acyl_smarts = "P(OC[C](O)COC(=O)[#6])"
    acyl_frag = Chem.MolFromSmarts(acyl_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl_frag, useChirality=False)
    if len(acyl_matches) != 1:
        if len(acyl_matches) == 0:
            return False, "No acyl ester substituent at the expected (sn-1) hydroxyl position found"
        else:
            return False, f"Found {len(acyl_matches)} acyl ester substituents; expected exactly 1 on the sn-1 hydroxyl"
    
    # If all checks pass, classify as correct.
    return True, "Contains phosphoserine head group with an acyl substituent at the 1-hydroxy position"


# Quick testing (only executed when running this script directly):
if __name__ == "__main__":
    test_smiles = [
        # Expected true:
        "P(OCC(O)COC(=O)CCCCCCC/C=C/CCCCCCCC)(OCC(N)C(O)=O)(O)=O",  # 1-Oleoylglycerophosphoserine
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O",  # PS(16:0/0:0)
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O"       # PS(12:0/0:0)
    ]
    for sm in test_smiles:
        result, reason = is_1_acyl_sn_glycero_3_phosphoserine(sm)
        print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")