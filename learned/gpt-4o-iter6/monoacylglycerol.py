"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    A monoacylglycerol consists of a glycerol backbone in which one hydroxy group
    is esterified with a fatty acid, and the other two
    hydroxy groups can be either free or substituted.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a glycerol backbone which is C-C-C with hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for an ester linkage (C(=O)O) 
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # We assume that a long hydrocarbon chain is an acyl group
    acyl_group_pattern = Chem.MolFromSmarts("C(=O)[C;!R]")
    if not mol.HasSubstructMatch(acyl_group_pattern):
        return False, "No appropriate acyl group found attached to ester linkage"

    return True, "Contains glycerol backbone with one acyl group and two hydroxyl groups"

# Test example
# smiles = "O(C(=O)CCCCCCC/C=C\CCCCCCCC)CC(O)CO" # Example of a monoacylglycerol
# print(is_monoacylglycerol(smiles))