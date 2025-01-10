"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol contains a myo-inositol ring, a glycerophosphate group, and two fatty acid chains.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for myo-inositol ring structure
    inositol_pattern = Chem.MolFromSmarts("C1(C(O)C(O)C(O)C(O)C(O)O)C1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol ring structure found"

    # Check for phosphate connected to inositol and glycerol
    phosphate_glycerol_pattern = Chem.MolFromSmarts("COP(O)(O)OC1CCOC1(C(O)O)O")
    if not mol.HasSubstructMatch(phosphate_glycerol_pattern):
        return False, "No correct phosphate linkage found between inositol and glycerol backbone"
    
    # Validate for two ester linkages representing fatty acid chains
    ester_pattern = Chem.MolFromSmarts("COC(=O)[#6]")
    if mol.GetSubstructMatches(ester_pattern) < 2:
        return False, f"Expected 2 fatty acid chains, found {mol.GetSubstructMatches(ester_pattern)}"

    return True, "Contains myo-inositol ring with glycerophosphate and two esterified fatty acid chains"