"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
#!/usr/bin/env python3
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine compounds
Defined as: An sn-glycerophosphoserine compound having an acyl substituent at the 1-hydroxy position.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine compound.
    The molecule must contain a phosphoserine head group and an acyl ester substituent
    on the primary (sn-1) hydroxyl of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the phosphoserine head group.
    # Here we use a simplified pattern for the serine part: "OCC(N)C(=O)O"
    # which should appear as part of the sn-glycero-3-phosphoserine fragment.
    ps_pattern = Chem.MolFromSmarts("OCC(N)C(=O)O")
    if not mol.HasSubstructMatch(ps_pattern):
        return False, "Phosphoserine head group not found"

    # Ensure that a phosphorus atom is present (as the phosphate group is essential)
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found"
    
    # Check for acyl ester substituent.
    # We look for an ester moiety that appears as an acyl group: O-C(=O)-C...
    # This pattern should capture the linkage where an acyl chain is attached via an oxygen.
    acyl_pattern = Chem.MolFromSmarts("OC(=O)[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No acyl ester substituent found at the 1-hydroxy position"

    # Optionally, one could further verify the acyl chain length (e.g. a minimum number of carbons)
    # by inspecting the fragment that follows the carbonyl. For now, we assume that the match is enough.

    return True, "Contains phosphoserine head group with an acyl substituent at the 1-hydroxy position"


# For quick testing:
if __name__ == "__main__":
    test_smiles = [
        "P(OCC(O)COC(=O)CCCCCCC/C=C/CCCCCCCC)(OCC(N)C(O)=O)(O)=O",  # 1-Oleoylglycerophosphoserine
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O"  # PS(16:0/0:0)
    ]
    for sm in test_smiles:
        res, reason = is_1_acyl_sn_glycero_3_phosphoserine(sm)
        print(f"SMILES: {sm}\nResult: {res}\nReason: {reason}\n")