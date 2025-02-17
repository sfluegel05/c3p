"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
"""
Classifies: Phosphatidylcholine (A glycerophosphocholine bearing two acyl substituents at positions 1 and 2)
"""
from rdkit import Chem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    Phosphatidylcholine is defined as a glycero-3-phosphocholine that bears two acyl substituents 
    (esterified fatty acids) at positions 1 and 2.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the phosphocholine head group.
    #    We look for a substructure matching the typical phosphate group attached to choline.
    head_group_smarts = "[OP](=O)([O-])OCC[N+](C)(C)C"
    head_group = Chem.MolFromSmarts(head_group_smarts)
    if head_group is None:
        return False, "Error defining phosphocholine head group SMARTS"
    
    if not mol.HasSubstructMatch(head_group):
        return False, "Phosphocholine head group not found"
    
    # 2. Check for the two acyl ester groups.
    #    We look for ester bonds in the form “OC(=O)”.
    #    However, we must exclude the ester that binds the phosphate head group.
    ester_smarts = "OC(=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Error defining ester SMARTS"
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    acyl_ester_count = 0
    # Loop through all matches found.
    for match in ester_matches:
        # match[0] is the oxygen in "OC(=O)".
        o_atom = mol.GetAtomWithIdx(match[0])
        # Check if this oxygen is bound to a phosphorus. If yes, it is likely part of the phosphate head.
        attached_to_p = any(neighbor.GetAtomicNum() == 15 for neighbor in o_atom.GetNeighbors())
        if not attached_to_p:
            acyl_ester_count += 1

    if acyl_ester_count != 2:
        return False, f"Found {acyl_ester_count} acyl ester group(s), but exactly 2 are required for phosphatidylcholine"
    
    return True, "Contains phosphocholine head group with 2 acyl ester groups attached to a glycero-3-phosphocholine scaffold"
    
# Example usage (uncomment for testing)
# test_smiles = "C([C@@](COC(CCCCCCCCCCCCCCCCCCCCCCC)=O)(OC(CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)[H])OP([O-])(=O)OCC[N+](C)(C)C"
# result, reason = is_phosphatidylcholine(test_smiles)
# print(result, reason)