"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
#!/usr/bin/env python3
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine compounds
Defined as: An sn-glycerophosphoserine compound having an acyl substituent 
at the 1-hydroxy position.
The algorithm requires that the molecule contains a phosphorus atom and that,
from this phosphorus, exactly one branch matches a serine head group while exactly one branch 
matches the acyl ester substructure expected at the sn-1 (primary) hydroxyl position.
"""

from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines whether a molecule is a 1-acyl-sn-glycero-3-phosphoserine compound.
    Our algorithm checks for:
      (1) Presence of a phosphorus atom.
      (2) A phosphoserine head group defined by a branch from phosphorus matching:
             OC[C](N)C(O)=O
          (this corresponds to the serine moiety).
      (3) An acyl ester substituent on the sn-1 hydroxyl defined by a branch from phosphorus matching:
             OC[C](O)COC(=O)[#6]
          (this pattern forces an ester bond attached to a carbonyl which then leads to a carbon chain).
    Both SMARTS are anchored on phosphorus so that they only count if they are directly 
    attached as substituents.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple; True and an explanation if the molecule belongs to the class,
                     otherwise False with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (1) Check that a phosphorus atom is present.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphorus atom found"

    # (2) Define SMARTS patterns anchored on phosphorus.
    # The serine branch is expected to be: P(OC[C](N)C(O)=O)
    serine_smarts = "[#15](OC[C](N)C(O)=O)"
    serine_query = Chem.MolFromSmarts(serine_smarts)
    # The acyl branch is expected to be: P(OC[C](O)COC(=O)[#6])
    acyl_smarts = "[#15](OC[C](O)COC(=O)[#6])"
    acyl_query = Chem.MolFromSmarts(acyl_smarts)

    # (3) Check that a phosphoserine head group is present.
    serine_matches = mol.GetSubstructMatches(serine_query, useChirality=False)
    if not serine_matches:
        return False, "Phosphoserine head group not found"
    
    # (4) Check for exactly one acyl ester substituent (for the sn-1 hydroxyl).
    acyl_matches = mol.GetSubstructMatches(acyl_query, useChirality=False)
    if len(acyl_matches) != 1:
        if len(acyl_matches) == 0:
            return False, "No acyl ester substituent at the 1-hydroxy position found"
        else:
            return False, f"Found {len(acyl_matches)} acyl ester substituents; expected exactly 1 on the sn-1 hydroxyl"

    # If all conditions are met, classify as true.
    return True, "Contains phosphoserine head group with an acyl substituent at the 1-hydroxy position"


# Quick testing when running this script directly
if __name__ == "__main__":
    test_smiles = [
        # Expected true cases:
        "P(OCC(O)COC(=O)CCCCCCC/C=C/CCCCCCCC)(OCC(N)C(O)=O)(O)=O",  # 1-Oleoylglycerophosphoserine
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O",  # PS(16:0/0:0)
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O",       # PS(12:0/0:0)
    ]
    for sm in test_smiles:
        res, reason = is_1_acyl_sn_glycero_3_phosphoserine(sm)
        print(f"SMILES: {sm}\nResult: {res}\nReason: {reason}\n")