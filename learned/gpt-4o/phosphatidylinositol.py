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
    # Update pattern to consider flexible inositol configurations 
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C(O)1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol ring structure found"

    # Check for glycerophosphate linkage
    # Use a pattern that captures structural variations: Glycerol connected via phosphate to inositol
    glycerophosphate_pattern = Chem.MolFromSmarts("COP(O)(O)OC[C@@H](O)C1COC1")
    if not mol.HasSubstructMatch(glycerophosphate_pattern):
        return False, "No correct glycerophosphate linkage found"

    # Validate for two ester linkages
    # Usual ester linkage format: -OC(=O)-[C of fatty acid]
    ester_pattern = Chem.MolFromSmarts("COC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Expected 2 fatty acid chains, found {len(ester_matches)}"

    return True, "Contains myo-inositol ring with glycerophosphate and two esterified fatty acid chains"