"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.

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

    # Pattern to detect the glycerol backbone part of a phosphatidylglycerol
    # Generalized for possible stereo/isomer configurations and common modifications
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO[P](=O)(O)O)(CO)O")  # Search for glycerol-phosphate linkage
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)

    if not glycerol_matches:
        return False, "No glycerol backbone with phosphate linkage found"

    # Pattern for ester-linked long fatty acid chains
    ester_fatty_acid_pattern = Chem.MolFromSmarts("OC(=O)C")  # This captures ester linkage to a fatty acid chain
    ester_matches = mol.GetSubstructMatches(ester_fatty_acid_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester-linked fatty acids, found {len(ester_matches)}"

    # Ensure presence of at least one phosphatidyl group
    phosphatidyl_pattern = Chem.MolFromSmarts("COP(O)(=O)OC")
    phosphatidyl_matches = mol.GetSubstructMatches(phosphatidyl_pattern)
    if not phosphatidyl_matches:
        return False, "Missing phosphatidyl group"

    return True, "Contains glycerol backbone with phosphatidyl group and ester-linked fatty acid chains"