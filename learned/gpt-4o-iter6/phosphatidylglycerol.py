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
    
    # Updated relaxed pattern for the glycerol phosphate backbone
    glycerol_phosphate_pattern = Chem.MolFromSmarts("OCC(O)COP(=O)(O)OC")  # Allowing for variations in stereochemistry
    
    # Check for glycerol phosphate backbone
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol phosphate backbone found"
    
    # Ensure the presence of ester-linked fatty acids (minimum 2 ester groups expected)
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester-linked fatty acids, found {len(ester_matches)}"
    
    # Check for specific non-target groups (not required)
    choline_pattern = Chem.MolFromSmarts("N(C)(C)C")
    ethanolamine_pattern = Chem.MolFromSmarts("OCCN")
    
    # We actually want to just make sure no typical phospholipid headgroups are included.
    if mol.HasSubstructMatch(choline_pattern) or mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "Contains non-phosphatidylglycerol headgroup (choline or ethanolamine)"
    
    return True, "Contains glycerol phosphate backbone with phosphatidyl group and ester-linked fatty acid chains"