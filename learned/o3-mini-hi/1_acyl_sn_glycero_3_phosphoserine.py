"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
#!/usr/bin/env python3
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine compounds
Defined as: An sn-glycerophosphoserine compound having an acyl substituent at the 1-hydroxy position.
The algorithm requires that the molecule contain (1) a phosphoserine head group,
(2) a phosphorus atom, and (3) exactly one glycerol branch from phosphorus that
matches the expected pattern for an acyl ester substituent at the sn-1 hydroxyl.
"""

from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine compound.
    The molecule must contain a phosphoserine head group (as evidenced by a serine fragment
    in combination with a phosphate group) and an acyl ester substituent on the primary (sn-1)
    hydroxyl of the glycerol backbone. To reduce false positives, the acyl substituent is only
    accepted if it is attached via the expected glycerol branch from the phosphorus.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule belongs to the class, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check that the phosphoserine head group is present.
    # We require a serine fragment attached through an oxygen: either "OCC(N)C(O)=O" 
    # or a variant with an extra carbon chiral marker (to cover common notations).
    serine_smarts = "OC[C](N)C(O)=O"
    serine_frag = Chem.MolFromSmarts(serine_smarts)
    if not mol.HasSubstructMatch(serine_frag):
        return False, "Phosphoserine head group not found"
    
    # 2. The molecule must contain a phosphorus atom.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found"
    
    # 3. Check for the acyl ester substituent at the sn-1 hydroxyl.
    # Here we require that one substituent from phosphorus exactly matches the expected branch:
    # The expected branch for the acyl substituent is a glycerol fragment with the acyl part:
    #   P(OC[C](O)COC(=O)R)
    # This SMARTS means: from phosphorus, an oxygen which is linked to a carbon (could be CH2 or CH)
    # that bears an OH, then a CH2 and then an oxygen that is esterified (OC(=O)) to an aliphatic (atom #6).
    acyl_smarts = "P(OC[C](O)COC(=O)[#6])"
    acyl_frag = Chem.MolFromSmarts(acyl_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl_frag)
    if len(acyl_matches) != 1:
        if len(acyl_matches) == 0:
            return False, "No acyl ester substituent at the expected (sn-1) hydroxyl position found"
        else:
            return False, f"Found {len(acyl_matches)} acyl ester substituents; expected exactly 1 on the sn-1 hydroxyl"
    
    return True, "Contains phosphoserine head group with an acyl substituent at the 1-hydroxy position"


# For quick testing:
if __name__ == "__main__":
    test_smiles = [
        "P(OCC(O)COC(=O)CCCCCCC/C=C/CCCCCCCC)(OCC(N)C(O)=O)(O)=O",  # 1-Oleoylglycerophosphoserine
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O"  # PS(16:0/0:0)
    ]
    for sm in test_smiles:
        result, reason = is_1_acyl_sn_glycero_3_phosphoserine(sm)
        print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")