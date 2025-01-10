"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    A monoacylglycerol has a glycerol backbone with one acyl group esterified at one of the hydroxyl positions,
    while the other two hydroxyl groups might remain free or chemically substituted.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string" 

    # Identify glycerol backbone - C1(CO)C(O)C1O pattern
    glycerol_pattern = Chem.MolFromSmarts("C(CO)(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for exactly one ester linkage - C(=O)O pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Validate long chain presence - typical acyl group, avoiding rings
    acyl_group_pattern = Chem.MolFromSmarts("C(=O)C-[C;!R]")
    if not mol.HasSubstructMatch(acyl_group_pattern):
        return False, "No appropriate acyl group found attached to ester linkage"

    return True, "Contains glycerol backbone with one acyl group and two hydroxyl groups"

# Example Test
# smiles = "O(C(=O)CCCCCCC/C=C\CCCCCCCC)CC(O)CO" # Example of a monacylglycerol
# print(is_monoacylglycerol(smiles))