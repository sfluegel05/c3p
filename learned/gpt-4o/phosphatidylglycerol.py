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
    
    # Pattern for glycerol backbone (C-C-C with 3 oxygens)
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](O)C(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Pattern for phosphatidyl group (phosphate with oxygen linkage)
    phosphate_pattern = Chem.MolFromSmarts("OC[P](=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing phosphatidyl group"
    
    # Pattern for ester linkages to fatty acids (O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester links, need at least 2"
    
    return True, "Contains glycerol backbone with attached phosphatidyl group and fatty acid ester links"