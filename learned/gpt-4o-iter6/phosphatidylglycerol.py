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
    
    # Updated pattern for the glycerol backbone with phosphate linkage
    # Here we enhance the pattern to ensure flexibility for different stereochemistry and connectivity
    glycerol_phosphate_pattern = Chem.MolFromSmarts("C([C@@H](COP(O)(=O)O)O)COP(=O)(O)OC")
    
    # Check for glycerol backbone with phosphate linkage
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with correct phosphate linkage found"
    
    # Ensure the ester linkage for the fatty acids (minimum 2 ester groups expected)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester-linked fatty acids, found {len(ester_matches)}"
    
    # Check to ensure specific non-target groups like choline/ethanolamine are absent
    choline_pattern = Chem.MolFromSmarts("N(C)(C)C")
    ethanolamine_pattern = Chem.MolFromSmarts("OCCN")
    if mol.HasSubstructMatch(choline_pattern):
        return False, "Contains choline group indicating not a phosphatidylglycerol"
    if mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "Contains ethanolamine group indicating not a phosphatidylglycerol"

    return True, "Contains glycerol backbone with phosphatidyl group and ester-linked fatty acid chains"