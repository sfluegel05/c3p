"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol is characterized by a glycerol backbone, phosphatidyl group, and fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Broadened pattern for glycerol backbone
    # Adjust to capture configurations and additional glycerol features
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)C(O)P")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Improved pattern for phosphatidyl group
    # This pattern encompasses variations with phosphate linkage
    phosphatidyl_pattern = Chem.MolFromSmarts("P(=O)(O)OC[C@@H](O)CO")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "Missing phosphatidyl group"
    
    # Check for at least two ester groups, which are typical in phosphatidylglycerol
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester links, need at least 2"
    
    return True, "Contains glycerol backbone with phosphatidyl group and fatty acid ester links"