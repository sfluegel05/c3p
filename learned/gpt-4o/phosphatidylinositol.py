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

    # Updated pattern for myo-inositol ring to allow for stereochemical variation
    inositol_pattern = Chem.MolFromSmarts("C1(CO)C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol ring structure found"

    # Check for glycerophosphate linkage with correct stereochemistry
    glycerophosphate_pattern = Chem.MolFromSmarts("COP(O)(O)OC[C@H]()[O][C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glycerophosphate_pattern):
        return False, "No correct glycerophosphate linkage found"

    # Validate for two ester linkages
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")  # Extend to match esterified fatty acid carbons
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Expected 2 fatty acid chains, found {len(ester_matches)}"

    return True, "Contains myo-inositol ring with glycerophosphate and two esterified fatty acid chains"