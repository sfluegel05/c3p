"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
"""
Classifies: Phosphatidylcholine
A phosphatidylcholine is defined as a glycero-3-phosphocholine bearing two acyl substituents 
at positions 1 and 2 (i.e. a phosphocholine head group attached to a glycerol esterified with 2 fatty acids).
"""
from rdkit import Chem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    Phosphatidylcholine is defined as a glycero-3-phosphocholine with two acyl ester groups
    (fatty acids) at positions 1 and 2 and a phosphocholine head group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a phosphatidylcholine, False otherwise.
        str: A description of the classification reasoning.
    """
    # Parse the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the phosphocholine head group by matching its characteristic substructure.
    # Changed SMARTS from "[OP](=O)([O-])OCC[N+](C)(C)C" to one that targets the phosphorus center.
    head_group_smarts = "P(=O)([O-])OCC[N+](C)(C)C"
    head_group = Chem.MolFromSmarts(head_group_smarts)
    if head_group is None:
        return False, "Error defining phosphocholine head group SMARTS"
    
    if not mol.HasSubstructMatch(head_group):
        return False, "Phosphocholine head group not found"
    
    # 2. Check for the two acyl ester groups.
    # We search for the ester bond fragment "OC(=O)". However, we must exclude the ester that connects 
    # the phosphate to the glycerol (i.e. the one in the phosphocholine head group).
    ester_smarts = "OC(=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Error defining ester SMARTS"
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    acyl_ester_count = 0
    # Loop through all ester matches.
    for match in ester_matches:
        # match[0] is the oxygen in "OC(=O)".
        o_atom = mol.GetAtomWithIdx(match[0])
        # Exclude if this oxygen is attached to a phosphorus (likely part of the head group).
        attached_to_p = any(neighbor.GetAtomicNum() == 15 for neighbor in o_atom.GetNeighbors())
        if not attached_to_p:
            acyl_ester_count += 1

    if acyl_ester_count != 2:
        return False, f"Found {acyl_ester_count} acyl ester group(s); exactly 2 are required for phosphatidylcholine"
    
    # If both criteria are met, it is a phosphatidylcholine.
    return True, "Contains a phosphocholine head group with 2 acyl ester groups on a glycero-3-phosphocholine scaffold"

# Example usage:
# test_smiles = "C([C@@](COC(CCCCCCCCCCCCCCCCCCCCCCC)=O)(OC(CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)[H])OP([O-])(=O)OCC[N+](C)(C)C"
# result, reason = is_phosphatidylcholine(test_smiles)
# print(result, reason)