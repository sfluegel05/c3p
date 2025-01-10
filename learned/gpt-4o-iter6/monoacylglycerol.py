"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    A monoacylglycerol has a glycerol backbone with one acyl group esterified at one of the hydroxyl positions,
    while the other two positions can have hydroxyl or other substituents.

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

    # Identify possible glycerol backbone arrangements - C(CO)(CO)O with flexible substitutions
    glycerol_pattern_1 = Chem.MolFromSmarts("C(CO)CO")  # e.g., R-C(CO)(O)
    glycerol_pattern_2 = Chem.MolFromSmarts("C(CO)OC")  # e.g., OC-C(O)(O)
    if not (mol.HasSubstructMatch(glycerol_pattern_1) or mol.HasSubstructMatch(glycerol_pattern_2)):
        return False, "No glycerol backbone found"

    # Check for exactly one ester linkage - C(=O)O pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Validate long acyl chain presence - typical acyl group, avoiding rings
    acyl_group_pattern = Chem.MolFromSmarts("C(=O)C-[C;!R]")
    if not mol.HasSubstructMatch(acyl_group_pattern):
        return False, "No appropriate acyl group found attached to ester linkage"

    return True, "Contains glycerol backbone with one acyl group and variable substituents"

# Example Test
# smiles = "O(C(=O)CCCCCCC/C=C\CCCCCCCC)CC(O)CO" # Example of a monacylglycerol
# print(is_monoacylglycerol(smiles))