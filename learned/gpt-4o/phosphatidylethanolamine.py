"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    Phosphatidylethanolamines have a glycerophospholipid structure with two hydrophobic chains,
    a glycerol backbone, a phosphate group, and an ethanolamine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for glycerol backbone pattern with generic chiral carbon connectivity
    glycerol_pattern = Chem.MolFromSmarts("[C@H](COP(O)OCCN)(OC([CX3](=O))[CX3](=O))") 
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Check presence of phosphate group connected: O=P(O)(O)
    phosphate_pattern = Chem.MolFromSmarts("O=P(O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group identified"
    
    # Look for ethanolamine group
    ethanolamine_pattern = Chem.MolFromSmarts("OCCN")
    if not mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "No ethanolamine group found"

    # Ensure two ester linkages
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Expected 2 ester linkages, found {len(ester_matches)}"

    return True, "Molecule classified as phosphatidylethanolamine"

# Test cases for phosphatidylethanolamine would be added here for validation